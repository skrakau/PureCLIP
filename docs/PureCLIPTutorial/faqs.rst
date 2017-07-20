
Trouble Shooting - FAQs
====================================


1. `Loading reference ... terminate called after throwing an instance of 'seqan::ParseError'`
 `what():  Unexpected character 'M' found.`

 The reference sequence is expected to have an alphabet of {'A', 'C', 'G', 'T', 'N'}.
 The error occurs when the character 'M' is found within the sequence, as it is the case for some reference sequences.

2. When analyzing CLIP-seq data where peaks of interest have high read start counts at individual positions (> 250) (assuming duplicates were removed properly and only uniquely mapping reads were used), please use the branch 'truncCount_unlimited' for the moment. The master version will be adjusted soon.


.. Note::
    We currently work on PureCLIP to improve its usability and robustness, and are therefore happy about any feedback. If you run into any problems, please do not hesitate to contact us, we would be happy to help.      
