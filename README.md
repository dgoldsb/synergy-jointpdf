# README #

Firstly, this module implements the notion of a "joint probability distribution" among stochastic variables, as well as many common operations for them. Secondly, it implements many information-theoretical quantities such as Shannon entropy, mutual information, information synergy, and information-based optimization procedures and robustness tests.

### QUICK START ###

* Basic PDF operations.

```python
from jointpdf import JointProbabilityMatrix

# randomly generated joint probability mass function p(A,B) of 2 discrete stochastic variables, each having 3 possible values
p_AB = JointProbabilityMatrix(2,3)

# obtain the marginal distribution p(A) by marginalization
p_A = p_AB[0]

# obtain the conditional probability distribution p(B|A)
# [can also directly use .conditional_distributions(...)]
p_B_given_A = p_AB - p_A

# create the joint pmf again from summing a marginal with a conditional, since p(a,b) = p(a)*p(b|a)
assert p_AB == p_A + p_B_given_A
```

* Basic information-theoretical measures

```python
# Shannon entropy
p_AB.entropy()
p_A.entropy()

# mutual information between I(A:B)
p_AB.mutual_information([0], [1])

# append a third variable C which is deterministically computed from A and B, i.e., such that I(A,B:C)=H(C)
p_AB.append_redundant_variables(1)
p_ABC = p_AB

# compute the information synergy that C contains about A and B (takes a while)
p_ABC.synergistic_information([2], [0,1])
```