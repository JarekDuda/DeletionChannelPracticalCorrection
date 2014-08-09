DeletionChannelPracticalCorrection
==================================

Practical correction and capacity estimation of deletion channel using correction trees 

While we understand well information channels like [binary symmetric channel](http://en.wikipedia.org/wiki/Binary_symmetric_channel) (BSC: fixed independent probability of bit flip) or [erasure channel](http://en.wikipedia.org/wiki/Binary_erasure_channel) (BEC: some values are marked as unknown), we don't even know the capacity of the simplest synchronization channel – [deletion channel](http://en.wikipedia.org/wiki/Deletion_channel) (BDC): for each bit there is an independent probability (p) of removing it from the stream (e.g. 1010 -> 100).

This implementation uses large state (64 bit) analogue of sequential decoding (forward only), what required replacing convolutional codes with specially designed code. The same encoding was previously used for BSC - implementation, discussion, analytically finding Pareto coefficients for BSC and BEC [can be found here](https://indect-project.eu/correction-trees/). To get some intuitions, [here is interactive simulator for BSC](http://demonstrations.wolfram.com/CorrectionTrees/).

We create potentially huge tree of possible corrections. Each node corresponds to decoded sequence of bytes up to this point. Its branches correspond to the choice of the next decoded byte and the number of bits used for this purpose: 8 without correction or a smaller number when we test scenarios with some applied deletion. In each node there is tested uniformly distributed checksum: R redundancy bits out of 8 bits (rate is (8-R)/8 for R=1 .. 7). These bits always agree for the proper correction and disagree with 1-2^-R probability for improper ones, allowing to statistically reduce the growth of the tree. After detecting an error, the most probable looking node is chosen for further expansion (heap). 

For deletion channel, different deletion patterns can lead to the same corrected sequence - such branches should be removed, what is achieved here by checking if given state hasn't been already spotted for given bit position.

The statistical behavior of such random tree is well described by Pareto coefficient (c), saying that increasing twice the limit of the number of nodes, reduces the probability of failure approximately 2^c times. We could improve performance by using more resources (memory and time) as long as c>0. For BSC and BEC this coefficient can be found analytically and c = 0 corresponds to channel capacity there - it is the boundary where the tree of possible corrections starts growing exponentially. 

It is argued that standard codebooks are not optimal for deletion channel, especially that they cannot work for p>1/2, while there is known (1-p)/9 universal lower rate bound ([survey paper](http://www.eecs.harvard.edu/~michaelm/TALKS/DelSurvey.pdf)). So roughtly extrapolated c=0 positions should be rather seen as lower bounds for low deletion probabilities here. For large deletion probability there are used codes with long sequences of the same value (0 or 1). However, in practical applications, deletions are rather low probable and appear alongside other types of damages like bit-flips, which would damage the block structure. In contrast, other types of errors can be easily added to the presented approach as just different types of branches with corresponding probabilities.

The tests for R=1, 4, 6, 7 (rate = 7/8, 1/2, 1/4, 1/8) were made for length 1000 byte encoded sequences (frames), with 5*10^7 node limit. 1000 frames were tested for each case. "damaged" is the number of improperly corrected frames out of 1000. "nodes" is the average number of created tree nodes per encoded byte - linear coefficient for time and memory cost (1 if error-free). Pareto coefficient (c) was estimated by linear fit to the central data ([1/3,2/3]). The last column contains roughtly extrapolated c=0 probability:

<table>
  <tr>
    <th>rate 7/8</th><th>p=0.001</th><th>0.003</th><th>0.005</th><th>0.007</th><th>0.009</th><th>0.011</th><th>0.012</th><th>0.013</th><th>0.014</th><th>~0.015</th>
  </tr>
  <tr>
    <th> c </th><th>8.85</th><th>1.54</th><th>0.65</th><th>0.346</th><th>0.245</th><th>0.138</th><th>0.114</th><th>0.071</th><th>0.030</th><th>0</th>
  </tr>
  <tr>
    <th> damaged </th><th>0</th><th>0</th><th>1</th><th>32</th><th>128</th><th>425</th><th>595</th><th>753</th><th>882</th><th>-</th>
  </tr>
  <tr>
    <th> nodes </th><th>1.23</th><th>10.5</th><th>235</th><th>2355</th><th>9602</th><th>25099</th><th>33823</th><th>40992</th><th>46054</th><th>-</th>
  </tr>
  <tr>
   <th></th><th></th><th></th><th></th><th></th><th></th><th></th><th></th><th></th><th></th><th></th>
  </tr>
  <tr>
    <th>rate 1/2</th><th>p=0.01</th><th>0.02</th><th>0.03</th><th>0.04</th><th>0.05</th><th>0.06</th><th>0.07</th><th>0.08</th><th>0.09</th><th>~0.1</th>
  </tr>
  <tr>
    <th> c </th><th>71</th><th>22</th><th>7.9</th><th>3.3</th><th>1.46</th><th>0.74</th><th>0.398</th><th>0.241</th><th></th><th>0</th>
  </tr>
  <tr>
    <th> damaged </th><th>0</th><th>0</th><th>0</th><th>0</th><th>0</th><th>1</th><th>22</th><th>170</th><th></th><th>-</th>
  </tr>
  <tr>
    <th> nodes </th><th>1.03</th><th>1.11</th><th>1.35</th><th>2.22</th><th>9.01</th><th>156</th><th>2093</th><th>12696</th><th></th><th></th>
  </tr>
  <tr>
    <th></th><th></th><th></th><th></th><th></th><th></th><th></th><th></th><th></th><th></th><th></th>
  </tr>
  <tr>
    <th>rate 1/4</th><th>p=0.1</th><th>0.12</th><th>0.13</th><th>0.14</th><th>0.15</th><th>0.16</th><th>0.17</th><th>0.18</th><th>0.19</th><th>~0.2</th>
  </tr>
  <tr>
    <th> c </th><th>7</th><th>2.9</th><th>2.1</th><th>1.4</th><th>0.89</th><th>0.66</th><th>0.46</th><th>0.287</th><th>0.153</th><th>0</th>
  </tr>
  <tr>
    <th> damaged </th><th>0</th><th>0</th><th>0</th><th>0</th><th>0</th><th>7</th><th>35</th><th>125</th><th>446</th><th>-</th>
  </tr>
  <tr>
    <th> nodes </th><th>1.54</th><th>2.80</th><th>4.88</th><th>11.7</th><th>84.7</th><th>746</th><th>2818</th><th>9996</th><th>27516</th><th>-</th>
  </tr>
  <tr>
   <th></th><th></th><th></th><th></th><th></th><th></th><th></th><th></th><th></th><th></th><th></th>
  </tr>
   <tr>
    <th>rate 1/8</th><th>p=0.1</th><th>0.12</th><th>0.14</th><th>0.16</th><th>0.18</th><th>0.20</th><th>0.22</th><th>0.24</th><th>0.26</th><th>?</th>
  </tr>
  <tr>
    <th> c </th><th>35.8</th><th>23.7</th><th>15.7</th><th>8.70</th><th>5.51</th><th>3.26</th><th>1.77</th><th>0.87</th>
    <th></th><th>0</th>
  </tr>
  <tr>
    <th> damaged </th><th>0</th><th>0</th><th>0</th><th>0</th><th>0</th><th>0</th><th>0</th><th>0</th><th></th><th>-</th>
  </tr>
  <tr>
    <th> nodes </th><th>1.10</th><th>1.19</th><th>1.20</th><th>1.40</th><th>1.82</th><th>2.92</th><th>9.40</th><th>97</th><th></th><th>-</th>
  </tr>
</table>

Estimated c=0 positions are close to known theoretical capacity bounds ([some recent article](http://arxiv.org/pdf/1211.2497v1.pdf). Comparing to other implementations, [here is some LDPC-based 2003 article](http://www.eecs.harvard.edu/~chaki/doc/code-long.pdf) which e.g. breaks for p~0.07-0.08 for rate 0.2333 code (page 8), while presented implementation still works above it for rate 1/2 - allowing to transmit more than twice more information through the same channel.

As this implementation uses the last state for final verification, sending this last state (protected) means that the rates should be reduced by a tiny factor, which decreases proportionally to frame length (Pareto coefficient does not depend on it). Without using this state, the last part of the message may remain damaged. Having this state, we can add bidirectional correction: simultaneously build tree in backward direction and finally merge both of them. This way we need two critical error concentrations to essentially stop the correction process, making that probability of failure is approximately squared (Pareto coefficient is doubled), what is confirmed by results for BSC. 

Another further work is trying different codewords (especially for large deletion probabilities) - the current coding was designed for bit-flips. Also, we can try to analytically find Pareto coefficients here – the difficulty in comparison to BSC is cutting branches (deletion patterns) corresponding to the same correction.

Jarek Duda, Kraków, 8/08/2014
