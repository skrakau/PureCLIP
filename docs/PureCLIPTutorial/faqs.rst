
Trouble Shooting - FAQs
====================================


1. `Loading reference ... terminate called after throwing an instance of 'seqan::ParseError'`
  `what():  Unexpected character 'M' found.` 

The reference sequence is expected to have an alphabet of {'A', 'C', 'G', 'T', 'N'}.
The error occurs when the character 'M' is found within the sequence, as it is the case for some reference sequences.


