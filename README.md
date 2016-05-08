# Unsupervised Learning of Mixture of von Mises–Fisher distribution

When fitting a data by a mixture model, the first problem we face is how to choose the number of mixture (compoents). This repository contains matlab code for unsupervised learning of mixture of <a href="https://www.wikiwand.com/en/Von_Mises%E2%80%93Fisher_distribution"> von Mises–Fisher distribution</a>. The best number of componets that fit the input data is selected based on <a href="https://www.wikiwand.com/en/Minimum_message_length">Minimum message length (MML) Criterion </a>.

# Description
The matlab codes are heavily commented (lines that need some explanation). Read the function help to know what goes in out and the function does. The two important files are ``` fitVMM_CEM.m ``` and ``` fitVMM_EM.m```. The first one fits the input data by performing the model selection (unsupervised learning). The second one need fits the data using the EM algorithm, hence need the number of mixture components as input to function.   

# Usage
 Run ```demo.m``` file. The ```demo-data.mat``` file contains a data sampled from mixture of three von Mises compoents. Below is the plot of this data 
<p align="center">
<img src="https://github.com/isrish/VMM/blob/master/transparent.png" height="260" width="240"/>
</p>

# TODO
