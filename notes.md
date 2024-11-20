## Data:

development and testing datasets:

- spa.covered.dev (root + features) --> spa.dev (answers produced by model)

- spa.covered.tst (root + features) --> spa.tst (answers produced by model)


training data: 

- spa.trn = training data with actual answers

    Each line consists of root ("lemma") + features + inflected form 


## nonneural.py:

base code for model (description on page 6 of https://aclanthology.org/2023.sigmorphon-1.13.pdf)

- Alignment: For each line in the training data,

-- aligns the root with the inflected form by finding minimum-cost edit path (substitutions higher cost than insertion/deletions)

-- splits each aligned pair into prefix + stem + suffix based on the central part that should match


- Inflection rules: given the aligned pairs, creates prefix and suffix changing rules based on how the prefix/suffix changed (suffix rules are also in format root+suffix -> root+modified suffix in addition to just suffix -> modified suffix)

- associate these rules with the relevant features 


- Generation: given lemma and features, applies rules associated with those features. No generalization.

-- Specifically, uses longest-matching suffix rule 




## diff file format

command: diff ./data/spa.dev ./data/spa.out > diff vN

format is 

start,end (of lines that were different)

< lines from .dev

---

> lines from .out