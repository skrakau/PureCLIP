
Trouble Shooting - FAQs
====================================


1. `Loading reference ... terminate called after throwing an instance of 'seqan::ParseError'`
 `what():  Unexpected character 'M' found.`

 The reference sequence is expected to have an alphabet of {'A', 'C', 'G', 'T', 'N'}.
 The error occurs when the character 'M' is found within the sequence, as it is the case for some reference sequences.


.. Note::
    We currently work on PureCLIP to improve its usability and robustness, and are therefore happy about any feedback. If you run into any problems, please do not hesitate to contact us, we would be happy to help.      
