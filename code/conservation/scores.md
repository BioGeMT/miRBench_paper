# Understanding Conservation Scores: phyloP and phastCons 100way

## What Are Conservation Scores?

Conservation scores are numerical values that measure how well DNA sequences have been preserved across species throughout evolution. These scores are calculated by:
1. Aligning genome sequences from multiple species
2. Comparing the rate of sequence changes
3. Calculating the probability that sequence changes are constrained by evolution

High conservation typically indicates:
* Functional importance of the sequence
* Selection pressure to maintain the sequence
* Essential biological roles

Low conservation or negative scores may indicate:
* Regions under positive selection
* Functionally neutral sequences
* Species-specific adaptations

## Overview of the 100way Scores

The "100way" refers to a multiple sequence alignment comparing the human genome (hg38) with 99 other vertebrate species. This alignment:
* Spans diverse vertebrate species
* Provides deep evolutionary context
* Enables sensitive detection of conservation
* Includes mammals, birds, fish, and reptiles

## phyloP (phylogenetic P-values)
**File**: `hg38.phyloP100way.bw`

### What is phyloP?
phyloP measures nucleotide-by-nucleotide evolutionary conservation by:
* Computing evolution rate at each base
* Comparing to a neutral model of evolution
* Testing for deviation from neutral drift
* Identifying both conservation and acceleration

### Score Calculation
* Based on likelihood ratio tests
* Considers phylogenetic tree structure
* Accounts for base substitution patterns
* Computes deviation from neutral evolution rate

### Score Interpretation
* **Scale**: -20 to +20 (typical range)
  * **Positive scores** (0 to +20)
    * Slower evolution than expected
    * Indicates purifying selection
    * Higher scores = stronger conservation
    * Common in functional elements
  
  * **Negative scores** (-20 to 0)
    * Faster evolution than expected
    * Indicates positive selection
    * Lower scores = more rapid change
    * May indicate adaptive elements

## phastCons (Phylogenetic Analysis with Space/Time Models Conservation)
**File**: `hg38.phastCons100way.bw`

### What is phastCons?
phastCons estimates conservation using a hidden Markov model (HMM) approach by:
* Modeling evolutionary processes
* Identifying conserved elements
* Computing conservation probabilities
* Considering regional patterns

### Score Calculation
* Uses two-state HMM
* Models conserved vs. non-conserved states
* Considers neighboring base influences
* Accounts for phylogenetic relationships
* Incorporates background mutation rates

### Score Interpretation
* **Scale**: 0 to 1 (probability)
  * **0**: No conservation detected
    * Evolves at neutral rate
    * No evidence of selective pressure
  
  * **1**: Complete conservation
    * Strong evolutionary constraint
    * High functional importance
  
  * Score ranges:
    * 0.8 - 1.0: Highly conserved (functional elements)
    * 0.5 - 0.8: Moderately conserved (possible function)
    * 0.0 - 0.5: Weakly conserved (likely neutral)

## Comparing phyloP vs phastCons

### Main Differences
1. **Resolution**
   * **phyloP**: Single-nucleotide resolution, independent sites
   * **phastCons**: Regional pattern detection, considers neighboring bases

2. **Methodology**
   * **phyloP**: Statistical test for rate deviation
   * **phastCons**: Probability model with state transitions

3. **Evolutionary Signal**
   * **phyloP**: Detects both conservation and acceleration
   * **phastCons**: Focuses on conservation detection

### Standard Thresholds
* **phyloP**: |2| typically used as significance threshold
* **phastCons**: 0.7 typically used as conservation threshold