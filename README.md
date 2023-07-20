# Multiplication-Free Biometric Recognition for Faster Processing under Encryption

## Description

This repository contains the implementation for generating the lookup tables of the pre-computed inner product (MFIP) and the pre-computed squared Euclidean distance (MFSED), freeing the IP and the SED from multiplication.
It also contains the proof-of-concept implementation of the integration of MFIP and MFSED with homomorphic encryption (HE) to perform biometric verification, which we call multiplication-free biometric recognition (MFBR).


Assuming normalized feature vectors, the MFIP (resp. MFSED) lookup table is parametrized by the feature vectors' dimension $d$, a feature quantization level $2^n$ where $n$ expresses the number of bits, and a cell quantization step $\Delta$, which we denote as $\text{MFIP}(d,n,\Delta)$ (resp.
$\text{MFSED}(d,n,\Delta)$).

## Remark

The use of the MFIP and MFSED lookup tables is not restricted to biometric recognition.
They can be used in any other application involving the computation of the IP or the SED of two normalized vectors, not necessarily feature vectors.

## Dependencies per folder

### MFIP-MFSED-Tables folder

This folder contains a Python 3.9 implementation that requires the following packages:

- [`NumPy`](https://numpy.org/)  
- [`SciPy`](https://scipy.org/)

### MFBR-BTP folder

This folder contains a C++ implementation that requires the following libraries:

- [`PALISADE version v1.11.5`](https://gitlab.com/palisade/palisade-release)
- [`OpenMP`](https://www.openmp.org/)

## Datasets

- Synthetic normalized feature vectors can be generated by executing `MFIP-MFSED-Tables/genSynVectors.py`
- Facial feature vectors of dimension 512 from [VGGFace2](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=8373813) dataset to extract facial feature vectors of dimension 512 using [ResNet-100](https://openaccess.thecvf.com/content_cvpr_2016/papers/He_Deep_Residual_Learning_CVPR_2016_paper.pdf) trained by two different losses: one trained with [ArcFace](https://openaccess.thecvf.com/content_CVPR_2019/papers/Deng_ArcFace_Additive_Angular_Margin_Loss_for_Deep_Face_Recognition_CVPR_2019_paper.pdf) and another one trained with [CosFace](https://openaccess.thecvf.com/content_cvpr_2018/papers/Wang_CosFace_Large_Margin_CVPR_2018_paper.pdf).

## Experiments per folder

### MFIP-MFSED-Tables folder

The following experiments generate the tables

1) Run `genBordersMFIPandMFSED.py` to generate the borders and the lookup tables
2) Run `genSynVectors.py` to generate synthetic samples
3) Run `compIPvsMFIP.py` to test IP vs MFIP and `compSEDvsMFSED.py` to test SED vs MFSED over synthetic samples
4) Run `testFacialFeatures.py` to test IP vs MFIP and SED vs MFSED over facial features vectors




### MFBR-BTP folder

The following experiments consider the HE-based BTPs using PALISADE BFVrns as an HE scheme: 
> * The MFBR BTPs that integrate MFIP and MFSED with HE
> * The HE-based baseline for the IP and SED 
> * Boddeti's BTP [[B18]](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=8698601) that performs the IP over normalized vectors, which we re-implemented with PALISADE.


Before launching the below experiments, execute the following commands

```
sudo mkdir build
cd build
sudo cmake ..
```


#### Experiments for measuring the runtime 

The file name of experiments measuring the runtime is `exp<BTP><MODE>.cpp` where `BTP = {MFIP, MFSED, IPBaseline, SEDBaseline, IPBoddeti}` for the HE-based BTPs that can be run in two modes: the clear-text comparison with the threshold `MODE = ClearComp` or the encrypted comparison with the threshold `MODE = EncComp` using [[BHP+21]](https://ieeexplore.ieee.org/abstract/document/9585508).

To run the above experiments, execute the following commands where `nBits = {128, 192, 256}` corresponds to the security levels.

```
cd build
sudo make exp<BTP><MODE>.cpp
./exp<BTP><MODE> nBits
```


#### Experiments for measuring the template size 


The following experiments measure the template size of MFBR and the baseline BTPs by storing them in binary files.


- The MFBR reference template size is `serialTemplateMFBR.cpp` 
- The `serialProbeAndTemplateBaseline.cpp`

To run the above experiments, execute the following commands where `nBits = {128, 192, 256}` corresponding the security levels.

```
cd build
sudo make <experiment.cpp>
./experiment nBits
```

- The MFBR probe template size is `serialProbeMFBR.cpp`

To run the above experiments, execute the following commands.

```
cd build
sudo make <experiment.cpp>
./experiment
```



### Cleaning data

To clean up the generated data at once, execute the following
```
sudo make clean-data
```

## References

[[ B18 ]](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=8698601)
[[BHP+21]](https://ieeexplore.ieee.org/abstract/document/9585508)



## Bibtex Citation

```
@inproceedings{bassit2022multiplication,
  title={Multiplication-Free Biometric Recognition for Faster Processing under Encryption},
  author={Bassit, Amina and Hahn, Florian and Veldhuis, Raymond and Peter, Andreas},
  booktitle={2022 IEEE International Joint Conference on Biometrics (IJCB)},
  pages={1--9},
  year={2022},
  organization={IEEE}
}
```
