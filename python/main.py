import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from scipy.sparse import coo_matrix
from scipy.ndimage import zoom as scipy_zoom
from Bio.PDB import PDBParser
import Protread
import time
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from pathlib import Path
import gc
import threading
from tkinterdnd2 import DND_FILES, TkinterDnD


def get_ca_coordinates(pdb_file):
    """Extracts X, Y, Z coordinates of Alpha Carbon (CA) atoms."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    
    coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    coords.append(residue['CA'].get_coord())
    
    return np.array(coords, dtype=np.float64)


class ProteinContactMapGUI:
    """Main GUI application for protein contact map analysis."""
    
    def __init__(self, root):
        self.root = root
        self.root.title("Protein Contact Map Analyzer")
        self.root.geometry("1200x800")
        
        # Variables
        self.pdb_file = tk.StringVar()
        self.threshold = tk.DoubleVar(value=8.0)
        self.max_pixels = tk.IntVar(value=8000)
        self.auto_display = tk.BooleanVar(value=True)
        self.save_on_generate = tk.BooleanVar(value=False)
        
        # State variables
        self.current_contact_map = None
        self.current_img = None
        self.current_info = {}
        self.processing = False
        
        self.setup_ui()
        
    def setup_ui(self):
        """Set up the user interface."""
        # Main container
        main_container = ttk.Frame(self.root, padding="10")
        main_container.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Configure grid weights
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_container.columnconfigure(1, weight=1)
        main_container.rowconfigure(2, weight=1)
        
        # Left panel - Controls
        control_frame = ttk.LabelFrame(main_container, text="Controls", padding="10")
        control_frame.grid(row=0, column=0, rowspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), padx=(0, 10))
        
        self.setup_control_panel(control_frame)
        
        # Right panel - Top: Info display
        info_frame = ttk.LabelFrame(main_container, text="Information", padding="10")
        info_frame.grid(row=0, column=1, sticky=(tk.W, tk.E, tk.N), pady=(0, 10))
        
        self.setup_info_panel(info_frame)
        
        # Right panel - Middle: Progress
        progress_frame = ttk.Frame(main_container)
        progress_frame.grid(row=1, column=1, sticky=(tk.W, tk.E), pady=(0, 10))
        
        self.progress_var = tk.StringVar(value="Ready")
        self.progress_label = ttk.Label(progress_frame, textvariable=self.progress_var)
        self.progress_label.pack(fill=tk.X)
        
        self.progress_bar = ttk.Progressbar(progress_frame, mode='indeterminate')
        self.progress_bar.pack(fill=tk.X, pady=(5, 0))
        
        # Right panel - Bottom: Visualization
        viz_frame = ttk.LabelFrame(main_container, text="Contact Map Visualization", padding="10")
        viz_frame.grid(row=2, column=1, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        self.setup_visualization_panel(viz_frame)
        
    def setup_control_panel(self, parent):
        """Set up the control panel with file selection and parameters."""
        row = 0
        
        # File selection section
        file_section = ttk.LabelFrame(parent, text="PDB File", padding="10")
        file_section.grid(row=row, column=0, sticky=(tk.W, tk.E), pady=(0, 10))
        row += 1
        
        # Drag and drop area
        drop_frame = tk.Frame(file_section, bg='#e0e0e0', relief=tk.SUNKEN, bd=2, height=100)
        drop_frame.pack(fill=tk.X, pady=(0, 10))
        drop_frame.pack_propagate(False)
        
        drop_label = tk.Label(
            drop_frame, 
            text="📁 Drop PDB file here\nor click Browse button below",
            bg='#e0e0e0',
            fg='#666666',
            font=('Arial', 10)
        )
        drop_label.place(relx=0.5, rely=0.5, anchor=tk.CENTER)
        
        # Register drop target
        drop_frame.drop_target_register(DND_FILES)
        drop_frame.dnd_bind('<<Drop>>', self.on_drop)
        
        # Browse button
        browse_btn = ttk.Button(file_section, text="Browse for PDB File...", command=self.browse_file)
        browse_btn.pack(fill=tk.X, pady=(0, 5))
        
        # Selected file display
        file_display_frame = ttk.Frame(file_section)
        file_display_frame.pack(fill=tk.X)
        
        ttk.Label(file_display_frame, text="Selected:").pack(side=tk.LEFT)
        self.file_label = ttk.Label(file_display_frame, textvariable=self.pdb_file, 
                                     foreground='blue', font=('Arial', 9, 'italic'))
        self.file_label.pack(side=tk.LEFT, padx=(5, 0))
        
        # Parameters section
        param_section = ttk.LabelFrame(parent, text="Parameters", padding="10")
        param_section.grid(row=row, column=0, sticky=(tk.W, tk.E), pady=(0, 10))
        row += 1
        
        # Threshold
        threshold_frame = ttk.Frame(param_section)
        threshold_frame.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Label(threshold_frame, text="Distance Threshold (Å):").pack(anchor=tk.W)
        threshold_spinbox = ttk.Spinbox(
            threshold_frame,
            from_=1.0,
            to=50.0,
            increment=0.5,
            textvariable=self.threshold,
            width=10
        )
        threshold_spinbox.pack(anchor=tk.W, pady=(5, 0))
        ttk.Label(threshold_frame, text="(Typical: 8.0 Å)", font=('Arial', 8), 
                 foreground='gray').pack(anchor=tk.W)
        
        # Max pixels for saved images
        max_pixels_frame = ttk.Frame(param_section)
        max_pixels_frame.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Label(max_pixels_frame, text="Max Pixels (for saving):").pack(anchor=tk.W)
        max_pixels_spinbox = ttk.Spinbox(
            max_pixels_frame,
            from_=500,
            to=20000,
            increment=1000,
            textvariable=self.max_pixels,
            width=10
        )
        max_pixels_spinbox.pack(anchor=tk.W, pady=(5, 0))
        ttk.Label(max_pixels_frame, text="(Larger = better quality, bigger file)", 
                 font=('Arial', 8), foreground='gray').pack(anchor=tk.W)
        
        # Options section
        options_section = ttk.LabelFrame(parent, text="Options", padding="10")
        options_section.grid(row=row, column=0, sticky=(tk.W, tk.E), pady=(0, 10))
        row += 1
        
        ttk.Checkbutton(
            options_section,
            text="Auto-display contact map",
            variable=self.auto_display
        ).pack(anchor=tk.W, pady=2)
        
        ttk.Checkbutton(
            options_section,
            text="Save image after generation",
            variable=self.save_on_generate
        ).pack(anchor=tk.W, pady=2)
        
        # Action buttons
        button_section = ttk.Frame(parent)
        button_section.grid(row=row, column=0, sticky=(tk.W, tk.E), pady=(0, 10))
        row += 1
        
        self.generate_btn = ttk.Button(
            button_section,
            text="Generate Contact Map",
            command=self.generate_contact_map,
            style='Accent.TButton'
        )
        self.generate_btn.pack(fill=tk.X, pady=(0, 5))
        
        self.save_btn = ttk.Button(
            button_section,
            text="Save Current Map...",
            command=self.save_contact_map,
            state=tk.DISABLED
        )
        self.save_btn.pack(fill=tk.X, pady=(0, 5))
        
        ttk.Button(
            button_section,
            text="Clear",
            command=self.clear_visualization
        ).pack(fill=tk.X)
        
        # Help section
        help_section = ttk.LabelFrame(parent, text="Help", padding="10")
        help_section.grid(row=row, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        help_text = (
            "How to use:\n\n"
            "1. Drop a PDB file or click Browse\n"
            "2. Adjust parameters if needed\n"
            "3. Click Generate Contact Map\n"
            "4. Use mouse to zoom/pan the map\n\n"
            "Contact Map Info:\n"
            "• White pixels = contacts\n"
            "• Black pixels = no contact\n"
            "• Diagonal = residue pairs\n"
            "• Symmetric around diagonal"
        )
        help_label = tk.Text(help_section, wrap=tk.WORD, height=15, width=30)
        help_label.insert('1.0', help_text)
        help_label.config(state=tk.DISABLED, font=('Arial', 9))
        help_label.pack(fill=tk.BOTH, expand=True)
        
    def setup_info_panel(self, parent):
        """Set up the information display panel."""
        self.info_text = tk.Text(parent, height=6, wrap=tk.WORD, font=('Courier', 9))
        self.info_text.pack(fill=tk.BOTH, expand=True)
        self.info_text.insert('1.0', "No contact map generated yet.\n\nLoad a PDB file to begin.")
        self.info_text.config(state=tk.DISABLED)
        
    def setup_visualization_panel(self, parent):
        """Set up the matplotlib visualization panel."""
        # Create matplotlib figure
        self.fig = Figure(figsize=(8, 7), dpi=100)
        self.ax = self.fig.add_subplot(111)
        
        # Initial empty plot
        self.ax.text(0.5, 0.5, 'No contact map loaded\n\nGenerate a contact map to visualize',
                    ha='center', va='center', fontsize=12, color='gray',
                    transform=self.ax.transAxes)
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        
        # Create canvas
        self.canvas = FigureCanvasTkAgg(self.fig, parent)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Add toolbar
        toolbar_frame = ttk.Frame(parent)
        toolbar_frame.pack(fill=tk.X)
        toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
        toolbar.update()
        
    def on_drop(self, event):
        """Handle file drop event."""
        # Get the dropped file path
        file_path = event.data
        
        # Clean up the path (remove curly braces if present)
        file_path = file_path.strip('{}')
        
        # Check if it's a PDB file
        if file_path.lower().endswith('.pdb'):
            self.pdb_file.set(Path(file_path).name)
            self.current_pdb_path = file_path
            self.update_info(f"Loaded file: {Path(file_path).name}")
        else:
            messagebox.showerror("Invalid File", "Please drop a PDB file (.pdb extension)")
    
    def browse_file(self):
        """Open file browser to select PDB file."""
        file_path = filedialog.askopenfilename(
            title="Select PDB File",
            filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")]
        )
        
        if file_path:
            self.pdb_file.set(Path(file_path).name)
            self.current_pdb_path = file_path
            self.update_info(f"Loaded file: {Path(file_path).name}")
    
    def update_info(self, message):
        """Update the info text panel."""
        self.info_text.config(state=tk.NORMAL)
        self.info_text.delete('1.0', tk.END)
        self.info_text.insert('1.0', message)
        self.info_text.config(state=tk.DISABLED)
    
    def generate_contact_map(self):
        """Generate contact map in a separate thread."""
        if not hasattr(self, 'current_pdb_path') or not self.current_pdb_path:
            messagebox.showwarning("No File", "Please select a PDB file first")
            return
        
        if self.processing:
            messagebox.showinfo("Processing", "Already processing a file. Please wait.")
            return
        
        # Disable generate button
        self.generate_btn.config(state=tk.DISABLED)
        self.processing = True
        
        # Start progress indicator
        self.progress_var.set("Processing...")
        self.progress_bar.start(10)
        
        # Run in separate thread to keep GUI responsive
        thread = threading.Thread(target=self._generate_contact_map_thread)
        thread.daemon = True
        thread.start()
    
    def _generate_contact_map_thread(self):
        """Worker thread for contact map generation."""
        try:
            pdb_path = self.current_pdb_path
            threshold = self.threshold.get()
            
            self.root.after(0, lambda: self.progress_var.set("Loading PDB file..."))
            
            # Load coordinates
            ca_coords = get_ca_coordinates(pdb_path)
            n = len(ca_coords)
            
            self.root.after(0, lambda: self.progress_var.set(f"Generating contact map for {n} residues..."))
            
            if n > 30000:
                self.root.after(0, lambda: messagebox.showinfo(
                    "Large Protein", 
                    f"Processing {n} residues. This may take several minutes..."
                ))
            
            # Generate contact map
            start_time = time.perf_counter()
            rows, cols, n = Protread.generate_contact_map_sparse(ca_coords, threshold)
            
            del ca_coords
            gc.collect()
            
            num_contacts = len(rows)
            total_pairs = (n * n - n) // 2
            contact_density = (num_contacts / total_pairs) * 100
            
            end_time = time.perf_counter()
            duration = end_time - start_time
            
            # Create image
            self.root.after(0, lambda: self.progress_var.set("Generating visualization..."))
            
            img = np.zeros((n, n), dtype=np.uint8)
            img[rows, cols] = 255
            img[cols, rows] = 255
            
            # Store results
            self.current_img = img
            self.current_info = {
                'pdb_name': Path(pdb_path).name,
                'threshold': threshold,
                'n_residues': n,
                'n_contacts': num_contacts,
                'contact_density': contact_density,
                'duration': duration
            }
            
            # Create sparse matrix
            data = np.ones(len(rows), dtype=np.int8)
            contact_map = coo_matrix((data, (rows, cols)), shape=(n, n), dtype=np.int8)
            self.current_contact_map = contact_map + contact_map.T
            
            # Update GUI
            self.root.after(0, self._update_visualization)
            
            # Save if requested
            if self.save_on_generate.get():
                self.root.after(0, self.save_contact_map)
            
        except Exception as e:
            self.root.after(0, lambda: messagebox.showerror("Error", f"Error processing PDB file:\n{str(e)}"))
            self.root.after(0, lambda: self.progress_var.set("Error occurred"))
        finally:
            self.root.after(0, self._processing_complete)
    
    def _update_visualization(self):
        """Update the visualization with the generated contact map."""
        if self.current_img is None:
            return
        
        info = self.current_info
        img = self.current_img
        
        # Update info panel
        info_text = (
            f"PDB File: {info['pdb_name']}\n"
            f"Distance Threshold: {info['threshold']} Å\n"
            f"Number of Residues: {info['n_residues']}\n"
            f"Contacts Found: {info['n_contacts']:,}\n"
            f"Contact Density: {info['contact_density']:.2f}%\n"
            f"Processing Time: {info['duration']:.2f} seconds"
        )
        self.update_info(info_text)
        
        # Update visualization
        if self.auto_display.get():
            self.ax.clear()
            self.ax.imshow(img, cmap='binary', origin='lower', interpolation='nearest')
            
            title = (f"{info['pdb_name']}\n"
                    f"Threshold: {info['threshold']} Å | Residues: {info['n_residues']} | "
                    f"Contacts: {info['n_contacts']:,}")
            self.ax.set_title(title, fontsize=10)
            self.ax.set_xlabel("Residue Index", fontsize=9)
            self.ax.set_ylabel("Residue Index", fontsize=9)
            
            self.fig.tight_layout()
            self.canvas.draw()
        
        # Enable save button
        self.save_btn.config(state=tk.NORMAL)
        
        self.progress_var.set("Complete!")
    
    def _processing_complete(self):
        """Clean up after processing."""
        self.progress_bar.stop()
        self.generate_btn.config(state=tk.NORMAL)
        self.processing = False
    
    def save_contact_map(self):
        """Save the current contact map to a file."""
        if self.current_img is None:
            messagebox.showwarning("No Map", "No contact map to save. Generate one first.")
            return
        
        file_path = filedialog.asksaveasfilename(
            title="Save Contact Map",
            defaultextension=".png",
            filetypes=[("PNG files", "*.png"), ("All files", "*.*")],
            initialfile=f"contact_map_{self.current_info['pdb_name']}.png"
        )
        
        if file_path:
            try:
                self.progress_var.set("Saving image...")
                self._save_contact_map_with_ui(file_path)
                messagebox.showinfo("Success", f"Contact map saved to:\n{file_path}")
                self.progress_var.set("Save complete")
            except Exception as e:
                messagebox.showerror("Error", f"Error saving file:\n{str(e)}")
                self.progress_var.set("Save failed")
    
    def _save_contact_map_with_ui(self, output_path):
        """Save contact map with proper UI elements."""
        img = self.current_img
        info = self.current_info
        n = img.shape[0]
        max_pixels = self.max_pixels.get()
        
        # Determine if downscaling is needed
        if n > max_pixels:
            scale_factor = max_pixels / n
            img_scaled = scipy_zoom(img.astype(float), scale_factor, order=1)
            img_scaled = (img_scaled > 127).astype(np.uint8) * 255
            display_n = img_scaled.shape[0]
        else:
            img_scaled = img
            display_n = n
        
        # Create figure
        dpi = 100
        fig, ax = plt.subplots(figsize=(12, 10), dpi=dpi)
        
        ax.imshow(img_scaled, cmap='binary', origin='lower', interpolation='nearest')
        
        title = (f"Protein Contact Map: {info['pdb_name']}\n"
                f"Threshold: {info['threshold']} Å | Residues: {n} | "
                f"Contacts: {info['n_contacts']:,} | Density: {info['contact_density']:.2f}%")
        ax.set_title(title, fontsize=11, pad=10)
        
        ax.set_xlabel("Residue Index", fontsize=10)
        ax.set_ylabel("Residue Index", fontsize=10)
        
        if n != display_n:
            tick_positions = ax.get_xticks()
            tick_labels = [f"{int(pos * n / display_n)}" for pos in tick_positions]
            ax.set_xticklabels(tick_labels)
            ax.set_yticklabels(tick_labels)
            
            note_text = f"Note: Image downscaled from {n}×{n} to {display_n}×{display_n} pixels for file size"
            fig.text(0.5, 0.02, note_text, ha='center', fontsize=8, style='italic', color='gray')
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight', format='png')
        plt.close(fig)
    
    def clear_visualization(self):
        """Clear the current visualization and reset."""
        self.current_img = None
        self.current_contact_map = None
        self.current_info = {}
        
        self.ax.clear()
        self.ax.text(0.5, 0.5, 'No contact map loaded\n\nGenerate a contact map to visualize',
                    ha='center', va='center', fontsize=12, color='gray',
                    transform=self.ax.transAxes)
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.canvas.draw()
        
        self.update_info("Visualization cleared.\n\nLoad a PDB file to begin.")
        self.save_btn.config(state=tk.DISABLED)
        self.progress_var.set("Ready")


def main():
    """Main entry point for the GUI application."""
    root = TkinterDnD.Tk()
    
    # Set up custom style
    style = ttk.Style()
    style.theme_use('clam')
    
    # Configure custom button style
    style.configure('Accent.TButton', background='#0078d4', foreground='white')
    
    app = ProteinContactMapGUI(root)
    
    # Center window on screen
    root.update_idletasks()
    width = root.winfo_width()
    height = root.winfo_height()
    x = (root.winfo_screenwidth() // 2) - (width // 2)
    y = (root.winfo_screenheight() // 2) - (height // 2)
    root.geometry(f'{width}x{height}+{x}+{y}')
    
    root.mainloop()


if __name__ == "__main__":
    main()