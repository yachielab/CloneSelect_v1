# QUEEN cloning script for a large catalog of plasmids
**Samuel King | Yachie Lab | Last updated: January 2023**

This is a QUEEN script for constructing plasmids by Digestion-Ligation, Gibson Assembly, and Golden Gate Assembly for the plasmids described in this study.

## Setting up QUEEN
- We recommend using Anaconda as a distribution and coding your QUEEN scripts in JupyterLab.
- This tutorial https://biocircuits.github.io/chapters/00_setting_up_python_computing_environment.html describes very clearly how to set up your environment with Anaconda.
- Write out the exact cloning plan for your desired plasmid step-by-step, including template plasmid names, primer names, enzyme names, etc., before starting your QUEEN script.
- If designing many plasmids, a spreadsheet is useful to carry the information. Then, you can code your functions to work with the spreadsheet items.
- The full tutorial for QUEEN can be found here: https://github.com/yachielab/QUEEN

After setting up Anaconda and getting JupyterLab going, start your script like this:
```python
%matplotlib agg
#!pip install python-queen  #install QUEEN here!
#!conda install graphviz    #graphviz is not automatically installed by Anaconda and you need to install it yourself here
import sys
from QUEEN.queen import *
set_namespace(globals())
from QUEEN import cutsite as cs
import re
import pandas as pd
import numpy as np
pd.set_option('display.max_colwidth', None)
```

## QUEEN wrapper functions
Here we provide the following new functions that build upon regular expression and QUEEN to improve QUEEN's functionality for a large catalog of plasmids:
### Sequence search functions:
* `find_primer_binding`: Provides a QUEEN feature list of a primer binding sequence on a template when given primer and template QUEEN objects.
* `primer_binding_object`: Provides a QUEEN feature object of a (18bp+) primer binding sequence on a template when given primer and template QUEEN objects.
* `primer_binding_object_short`: Provides a QUEEN feature object of a short (15bp+) primer binding sequence on a template when given primer and template QUEEN objects.
* `find_primer_binding_str`: Provides a string of a primer binding sequence on a template when given primer (specified as 'FW' or 'RV") and template QUEEN objects.
* `primer_binding_start_pos`: Provides the start index position of the template sequence of a primer binding on a searched template sequence.
* `primer_binding_end_pos`: Provides the end index position of the template sequence of a primer binding on a searched template sequence.
* `find_primer_binding_primer`: Provides a QUEEN feature list of a primer binding sequence on another primer when given two primer QUEEN objects.
* `find_primer_binding_primer_str`: Provides a string of a primer binding sequence on another primer when given two primer QUEEN objects.
* `find_primer_binding_primer_start`: Provides the start index position of the binding portion of a primer sequence bound on a searched primer.
* `find_primer_binding_primer_end`: Provides the end index position of the binding portion of a primer sequence bound on a searched primer.
* `find_similar_sequence`: Provides the string of the specified binding portion of a primer that contains a specified number of mismatches to a template/plasmid string.

### DNA modification functions:
* `generate_complement`: Provides a string of the complementary sequence to the provided string DNA (not the reverse complement, this is what flipdna does).
* `generate_dsdna`: Provides a double-stranded QUEEN object when given a single-stranded string DNA.
* `anneal_oligos`: Creates a double-stranded QUEEN object with or without overhangs on either side when given two single-stranded QUEEN objects.
    * __Helper functions__:
        * `five_prime_overhang`: Provides the 5' overhang of a ssDNA with a known annealing sequence to another ssDNA.
        * `three_prime_overhang`: Provides the 3' overhang of a ssDNA with a known annealing sequence to another ssDNA.
        * `primer_overhangs`: Provides both the 5' and 3' overhangs of a ssDNA with a known annealing sequence to another ssDNA.
* `gibson_assembly`: Provides the circular ligated product QUEEN object of 2 overlapping DNA fragments (akin to a Gibson Assembly reaction).
* `stitch_fragments`: Provides the linear ligated product QUEEN object of 2 overlapping DNA fragments.
    
### Wrapper cloning functions:
* `create_pcr_product`: Provides a PCR product QUEEN object when given primer and template QUEEN objects.
* `create_pcr_product_special`: Provides a PCR product QUEEN object when given a *short* Fw primer (15bp+), Rv primer, and template QUEEN objects (special case).
* `create_pcr_product_mismatches`: Provides a PCR product QUEEN object when given a Fw primer with mismatches, a Rv primer with no mismatches, and a template/plasmid QUEEN object.
* `double_digest_insert`: Provides a digested QUEEN object when given two restriction enzymes and the insert QUEEN object.
* `double_digest_backbone`: Provides a digested QUEEN object when given two restriction enzymes and the backbone QUEEN object.
* `typeIIS_digest_insert`: Provides a digested QUEEN object when given a Type-IIS restriction enzyme and the insert QUEEN object.
* `typeIIS_digest_backbone`: Provides a digested QUEEN object when given a Type-IIS restriction enzyme and the backbone QUEEN object. 
