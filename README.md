# wsp

The Wedding Seating Problem is a combinatorial optimization problem incorporating elements of bin packing with conflicts, bin packing with profits, and load balancing.


## The Wedding Seating Problem

In the Wedding Seating Problem (WSP), groups of guests of different sizes must be seated at tables of limited capacities. Some of these groups may or may not like each other, thus some relation is defined over each pair of them. Pairs of groups whose relation is *definitely apart* can never be seated at the same table. While not strictly necessary, pairs of groups whose relation is either *rather together* or *rather apart* should, if possible, be seated together or apart, respectively. Pairs which have no specific relation are *indifferent*. Tables should also be balanced as much as possible.


## Instructions

1. Generate instances with `/tools/instance_generator.py`. Specific instructions can be found in the source code.

2. Solve instances with the various models, which can be found in `/models`.


## References

Refer to [this paper](http://cerc-datascience.polymtl.ca/wp-content/uploads/2018/01/Technical-Report_DS4DM-2017-015.pdf) for more details concerning the models.