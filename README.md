# BalanSiNG: Fast and Scalable Balanced Signed Network Generator Following Real-world Properties

## Overview
How can we efficiently generate large-scale signed networks following real-world properties?
In signed networks, relations between nodes are represented as positive and negative edges.
These signed networks have spurred great attention in data mining  communities, and many researchers have been engaged in developing meaningful applications on signed networks.
However, our understanding of signed network generation following realistic properties was immature.
Many models have been proposed for generating unsigned networks, but those models do not consider the desirable formation of signed edges at all.

In this paper, we propose `BalanSiNG`, a novel, scalable, and fully parallelizable method for generating large-scale signed networks following realistic properties.
We first identify a self-similar balanced structure observed from a real signed network, and then simulate the self-similarity via Kronecker product.
Furthermore, we propose noise and weight techniques to produce practical signed networks.
Then, we show that `BalanSiNG` efficiently produces signed edges fully in parallel.
Our experimental results demonstrate that `BalanSiNG` efficiently generates realistic signed networks satisfying various desired properties compared to its competitors.

## Paper
BalanSiNG is described in the following paper:
* BalanSiNG: Fast and Scalable Balanced Signed Network Generator Following Real-world Properties
    - Submitted to KDD 2019

## Datasets
The datasets used in this paper are avariable at this repository. 

## Usage
We provide two version of implementations based on `scala` and `c++`.

### How to use scala code

### How to use c++ code
