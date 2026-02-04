# Protein Contact Map Analyzer

I developed this application to provide a high-performance graphical interface for generating and visualizing protein contact maps from PDB files. By combining a responsive Python-based GUI with a high-performance C++ backend, I have ensured that even large proteins can be processed and visualized efficiently.

---

## Useful For?

I designed this tool for researchers and bioinformaticians who need to:

* Identify structural patterns and tertiary interactions within protein structures.
* Visualize residue-residue proximity within a specified distance threshold.
* Process massive protein structures that typically cause performance lag in standard Python implementations.
* Generate high-resolution contact map images for publication or further analysis.

---

## Installation Steps

To get this project running on your local machine, I recommend the following steps:

### 1. Prerequisites

Ensure you have a C++ compiler supporting C++11 and OpenMP, and a Python 3.x environment.

### 2. Compile the C++ Backend

I use a custom C++ extension named `Protread` to handle the heavy computations. You must compile the `Protread.cpp` file into a shared library or Python module compatible with your system.

```bash
# Example using g++ (adjust for your specific Python binding method)
g++ -O3 -shared -std=c++11 -fPIC -fopenmp $(python3 -m pybind11 --includes) Protread.cpp -o Protread$(python3-config --extension-suffix)

```

### 3. Install Python Dependencies

I have utilized several libraries for the GUI and data handling. Install them via pip:

```bash
pip install numpy matplotlib scipy biopython tkinterdnd2

```

---

## Usage

Once you have compiled the backend and installed the dependencies, you can launch the application:

1. **Start the App:** Run `python main.py`.
2. **Load a PDB:** Drag and drop a `.pdb` file directly into the interface or use the **Browse** button.
3. **Set Parameters:** Adjust the **Distance Threshold (Å)** (default is 8.0) and the **Max Pixels** for output quality.
4. **Generate:** Click **Generate Contact Map**. I have implemented threading so the GUI remains responsive during calculation.
5. **Explore:** Use the Matplotlib toolbar to zoom in on specific regions of the contact map.
6. **Save:** Export the map as a PNG file using the **Save Current Map** button.

---

## Optimizations Used

I prioritized performance by offloading the computational bottlenecks to C++. Below are the specific optimizations I implemented:

### Spatial Hashing

Instead of a naive  comparison where every residue is checked against every other residue, I implemented a **Spatial Hash Grid**.

* I partition the 3D space into cells with a size equal to the distance threshold.
* I only check residues within the same or neighboring cells (27 cells total).
* This reduces the average time complexity significantly, especially for large, sparse structures.

### Multi-threading with OpenMP

I utilized **OpenMP** to parallelize the contact search.

* The workload is distributed across all available CPU cores using a dynamic schedule.
* I use thread-local storage for contact pairs to minimize synchronization overhead, only merging them into the global list in a critical section at the end.

### Sparse Data Handling

Since most residue pairs are not in contact, I treat the contact map as a sparse matrix.

* I use `scipy.sparse.coo_matrix` to handle the data in Python, which saves massive amounts of memory compared to a dense  array.
* I use `shrink_to_fit()` in C++ to ensure memory usage is kept at a minimum after the vectors are populated.

### Memory Management

I explicitly call `gc.collect()` and delete large intermediate coordinate arrays in the Python wrapper to prevent memory leaks when processing multiple files in a single session.

---

### Tests

This was able to handle creating a protein contact map for "4PTH.pdb" (no interactive map).

Distance Threshold: 8.0 Å

Number of Residues: 40000
Contacts Found: 51,851,718
Contact Density: 6.48%
Processing Time: 1.45 seconds

---

## License

I provide this project under the MIT License. Feel free to use and modify it for your research needs.
