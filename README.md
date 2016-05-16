# Unsupervised Learning of Mixture of von Mises–Fisher distribution

When fitting data by a mixture model, the first problem we face is how to choose the number of mixture (components). This repository contains matlab code for unsupervised learning for  mixture of <a href="https://www.wikiwand.com/en/Von_Mises%E2%80%93Fisher_distribution"> von Mises–Fisher distribution</a>. The number of mixture component that best fit the input data is selected based on <a href="https://www.wikiwand.com/en/Minimum_message_length">Minimum Message Length (MML) criterion </a>.

# Description
The matlab codes are heavily commented (just lines that need some kind of explanation). Read the function help to know what goes in and out and what the function does. The two important files are ``` fitVMM_CEM.m ``` and ``` fitVMM_EM.m```. The first one fits the input data by performing model selection (unsupervised learning). The second one need fits the data using the EM algorithm, hence need the number of mixture components as an input.   

# Usage
 The ```demo-data.mat``` file contains simulated data sampled from mixture of three von Mises compoents. Run ```demo.m``` file. Below is the plot of demo data result: 
<p align="center">
<img src="https://github.com/isrish/VMM/blob/master/transparent.png" height="280" width="260"/>
</p>

# References
The code is based on ideas described in the following two papers: 
``` 
@inproceedings{banerjee2005clustering,
  title={Clustering on the unit hypersphere using von Mises-Fisher distributions},
  author={Banerjee, Arindam and Dhillon, Inderjit S and Ghosh, Joydeep and Sra, Suvrit},
  booktitle={Journal of Machine Learning Research},
  pages={1345--1382},
  year={2005}
}

@article{kasarapu2015minimum,
  title={Minimum message length estimation of mixtures of multivariate Gaussian and von Mises-Fisher distributions},
  author={Kasarapu, Parthan and Allison, Lloyd},
  journal={Machine Learning},
  volume={100},
  number={2-3},
  pages={333--378},
  year={2015},
  publisher={Springer}
}
```

