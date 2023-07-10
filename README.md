# SCAD SVM Classification with Weighted Data

This repository contains an R implementation of SCAD SVM (Support Vector Machine) classification with weighted data. SCAD SVM is a variant of SVM that incorporates variable selection using the SCAD (Smoothly Clipped Absolute Deviation) penalty. The implementation in this repository extends the SCAD SVM algorithm to handle weighted data, allowing for improved handling of imbalanced datasets.

## Table of Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Parameters](#parameters)

## Introduction

The SCAD SVM classification algorithm is a powerful approach for binary classification, combining the strengths of SVM with variable selection using the SCAD penalty. In this implementation, we have enhanced the algorithm by incorporating weighted data. One notable application of this approach is demonstrated in the paper titled "On Sparse representation for Optimal Individualized Treatment Selection with Penalized Outcome Weighted Learning".

## Installation

To use the SCAD SVM classification code with weighted data, follow these steps:

1. Ensure you have R installed on your system.
2. Clone this repository or download the code files.
3. Open R and set the working directory to the location where you saved the code files.

## Usage

To use the SCAD SVM classification function with weighted data, follow these steps:

1. Prepare your training data in the form of a data matrix `x` and corresponding target values `y`.
2. Optionally, provide class weights if you want to account for the imbalance in the dataset.
3. Call the SCAD SVM classification function and provide the necessary parameters, such as `lambda1`, `x`, `y`, `sample.weights`, and others.
4. The function will perform the SCAD SVM classification on the weighted data and return the result, including the weights, bias, selected feature indices, and other relevant information.

Example usage:

```R
# Prepare your data
x <- ...
y <- ...
sample.weights <- ...

# Call the SCAD SVM classification function
result <- scadsvc(lambda1 = 0.01, x = x, y = y, sample.weights = w_hat)

# Access the results
weights <- result$w
bias <- result$b
selected_indices <- result$xind
# ...
```

## Parameters
The SCAD SVM classification function with weighted data accepts several parameters, including:
- lambda1: Tuning parameter in the SCAD function.
- x: Data matrix of size n-by-d representing the training data (n samples, d features).
- y: Column vector of target values or class labels.
- sample.weights: Optional parameter for assigning weights to different samples.
- Additional parameters for customization (e.g., tolerance, seed, maxIter, verbose, etc.)
Please adjust the parameters according to your specific use case.

