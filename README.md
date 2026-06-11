# TrianglesEdgeDetection

Fast detection of faint curved edges in noisy images using a triangle-based binary partitioning scheme.

This repository provides a MATLAB implementation of the method presented in:

> **Fast Detection of Curved Edges at Low SNR**  
> Nati Ofir, Meirav Galun, Boaz Nadler, Ronen Basri  
> *IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2016*  
> https://www.cv-foundation.org/openaccess/content_cvpr_2016/html/Ofir_Fast_Detection_of_CVPR_2016_paper.html

---

## 📌 Overview

Detecting faint edges in noisy images is a fundamental problem in computer vision. Standard gradient-based methods often fail when edges are weak or corrupted by noise.

This project implements an efficient algorithm for detecting **curved edges under low signal-to-noise ratio (SNR)** conditions.

The method is based on:

- Partitioning the image into **triangular regions**
- Evaluating a **statistical edge score** inside each region
- Using a **hierarchical binary partitioning scheme** to efficiently search for edges

This approach enables robust detection of weak, curved structures that are difficult to recover using classical methods.

---

## 🧠 Method Highlights

- Triangle-based image decomposition
- Multi-scale edge detection
- Statistical scoring of candidate edges
- Efficient hierarchical search
- Robust to noise and low contrast

---

## ⚙️ Installation

No installation is required.

Clone the repository:

git clone https://github.com/NatiOfir/TrianglesEdgeDetection.git

---
## Usage

I = im2double(imread('img.png'));

result = runIm(I);

---
## Cite

@inproceedings{Ofir2016FastEdges,
  title={Fast Detection of Curved Edges at Low SNR},
  author={Ofir, Nati and Galun, Meirav and Nadler, Boaz and Basri, Ronen},
  booktitle={CVPR},
  year={2016}
}
