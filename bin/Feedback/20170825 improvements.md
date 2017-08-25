# First of all

"Ik begrijp de equal chance initialization niet."

"Ik zou het aantal hypothesen minimaliseren of verdelen in prioriteiten.

En zou je bij elke hypothese een heel korte uitleg willen geven hoe je dat bereikt? 

De eerste 2 en de 4e lijken me sowieso prima. De pareto zou super interessant zijn, ik zou eigenlijk twijfelen over de haalbaarheid maar je had geloof ik al wat optimalisatie geprogrammeerd?"

# Questions to Rick

* Look at the 'proof' that the complexity profile is non-decreasing

# Improved 23-8-2017

* Made the distinction in the introduction between zero memory and maximal resilience and maximal memory and zero resilience a bit more nuanced; my intuition is that a rigid system still has no real memory, as in my opinion memory involves keeping information from the initial state over time, regardless of the initial system state
* Split paragraph in the general introduction, and expanded the latter paragraph to be clearer
* I improved the definition for complexity in my introduction, and went over all mentions to make sure I was not confusing terms (I did refer to a different type of complexity a few times, such as ecosystem complexity, which I expressed in a different way)
	* I added the term 'complicated' as a placeholder for this, to make the distinction
* Added the quantification of memory and resilience to my first step
* Explained better why a smaller system size is more convenient (**note**: I need to rewrite this section anyway after my experiments are done, so I won't do a complete rewrite now)
* Added that synergy is the quantification of complexity in our opinion in the introduction

# Improve later

* Expand on that no synergy measure yet satisfies all axioms
* Scan the literature review on the distinction between complexity and complicatedness
* In the rewrite of the introduction of IT in ecology, take a careful look where I talk about complicatedness and complexity...  The two become intertwined here, as we attempt to measure complicatedness with the entropy, and later more advanced IT tools
* Go over the final part of my introduction to make sure it is consistent with the experiments I did in the end
	* Give much stronger motivation for the model that I chose
	* Computational complexity is a fine reason
	* Availability of real data as well
	* Defend my assumptions (such as that the Boolean function is determinstic) here as well
* Remove the first ... second ... third ... finally ... construction from my introduction (final paragraph)
	* Third hypothesis is optional, remove if I don't actually get to it
	* In the end it boils down to "In general, what is at minimum needed to show that the trade off between the two competing fitness functions (memory and resilience) leads to higher synergy? Note that this question is not yet biological. Second question: do biological systems tend to have more synergy? Well for that we would have to take indeed some real-world network data, distill a local input-output model for each node, and see if the synergy in these nodes is higher than randomly expected (for some convenient choice of defining 'randomly')."
* Shoot down older synergy measures
	* "Other approaches can then still be mentioned, but one for one they should be 'shot down' for having some undesirable property or not being exactly what we want." as synergy is the ultimate measure with the taken complexity definition
* Dive again into the pros and cons
* Fill all empty citations
* Rewrite availabable GRN models
* Check if I am consistent in the symbols I use, and if I don't use symbols for different things in subsequent formulas
* If I step away from continuous, remove from methods