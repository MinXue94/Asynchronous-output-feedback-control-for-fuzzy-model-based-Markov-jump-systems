# Asynchronous-output-feedback-control-for-fuzzy-model-based-Markov-jump-systems
This project allows you to solve a convex optimization problem and, hence, acieving the controller gians and the optimal dissipative performance. 
### Prerequisites
For running the project you need to instal Matlab 2016b or newer versions and LMI Toolbox. 
## Description
### The system model
We use an example of a single-link robot arm system. The system matrices, transition probability matrix, conditional probability matrix, dissipative parameters, probability of packet dropouts are given.
### The asynchronous output feedback controller
The designed controller has been built as an output feedback controller based on a hidden Markov model and a compensation scheme. 
The control gain matrices are computed by solving the onvex optimization problem of Theorem 2. Meanwhile, the optimal dissipative value can be obtained.
## Running the tests
For running this project you must include all the folder in the Matlab path. Then, just run
``` 
new_lmi.m 
```
Then, you can plot the figures of the closed-loop system state responses under the control input.
```
dissipativity_fig.m
```
