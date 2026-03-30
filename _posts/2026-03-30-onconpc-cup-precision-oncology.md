---
layout: post
title: "OncoNPC: Machine Learning for Cancer of Unknown Primary in Precision Oncology"
date: 2026-03-30
permalink: /blog/onconpc-cup-precision-oncology/
published: true
categories: [paper-review]
tags:
  - AI
  - precision-medicine
  - cancer-genomics
  - oncology
  - machine-learning
description: "A paper review of OncoNPC, a machine-learning model that uses targeted tumor sequencing to predict the primary site of cancers of unknown primary and guide treatment decisions."
---

Cancer of unknown primary (CUP) is one of the most frustrating diagnoses in oncology. A patient presents with metastatic disease, pathology confirms malignancy, but the original tissue of origin cannot be confidently assigned. That uncertainty matters because modern oncology is organized around **primary site**, **molecular subtype**, and increasingly **targetable genomic context**. If you do not know what the cancer is, it becomes much harder to know how to treat it.

This is why the 2023 Nature Medicine paper on **OncoNPC** is so interesting. Rather than asking clinicians to wait for ever more stains, scans, and expert review, the paper asks a practical precision-oncology question:

> Can targeted tumor sequencing data itself provide enough information to infer the likely primary cancer type, and can that prediction improve treatment decisions?

Formally, the model learns a multiclass predictor

$$
p_\theta(y = k \mid x),
$$

where $x$ represents genomic and clinical features extracted from a tumor sample and $y$ is one of the candidate primary cancer types. In this setting, the goal is not to prove the true tissue of origin with certainty. The goal is to generate a clinically useful posterior over plausible origins.

According to the [Nature Medicine article](https://www.nature.com/articles/s41591-023-02482-6), OncoNPC was trained on targeted next-generation sequencing data from **36,445 tumors across 22 cancer types from three institutions**. On held-out samples, the model achieved a **weighted F1 score of 0.942** for high-confidence predictions, and when applied to **971 CUP tumors**, it produced high-confidence predictions in **41.2%** of cases. Most importantly, the paper reported that patients with CUP who received first palliative treatments concordant with the OncoNPC prediction had better outcomes, with **hazard ratio = 0.348**.

That last point is what pushes this paper beyond a technical classifier benchmark. This is not just "AI for classification." This is a paper about **clinical decision support in precision medicine**.

Explore the full [Nature Medicine article](https://www.nature.com/articles/s41591-023-02482-6) and the authors' [OncoNPC code repository](https://github.com/itmoon7/onconpc) for the original work and implementation details.

{% include figure.liquid loading="eager" path="assets/img/blog/onconpc/figure1-pipeline.svg" class="img-fluid rounded z-depth-1" zoomable=true %}

<div class="caption">
    Figure 1. OncoNPC pipeline overview. Targeted tumor sequencing and basic clinical variables are converted into a confidence-aware posterior over likely primary cancer types.
</div>

---

## Why CUP Is a Precision Medicine Problem

The phrase _precision medicine_ is often used loosely, but CUP is actually one of the clearest use cases for it. In CUP, the bottleneck is not simply prognosis. It is **actionable classification**.

Traditional work-up includes:

- radiology
- pathology and immunohistochemistry
- patient history
- molecular profiling

Even after all of that, a subset of tumors remain unclassified. The clinical consequence is severe: many therapies are approved or chosen in a site-specific way, and targeted therapy matching becomes less straightforward when the tissue of origin is unknown.

The OncoNPC paper frames this problem well. Rather than replacing pathology, the model acts as an additional source of evidence derived from targeted tumor sequencing. That matters because panel sequencing is already widely performed in oncology practice. A model that runs on data clinicians already collect has a much lower barrier to adoption than one requiring a completely new assay.

---

## The Core Idea Behind OncoNPC

OncoNPC uses targeted sequencing and basic clinical variables to predict one of 22 primary cancer types. The model is built on **XGBoost**, and the article's code-availability section states that the authors used Python packages including `xgboost` and `shap` for model development and interpretation.

At a high level, the input vector can be thought of as

$$
x = [m, c, s, a, z],
$$

where:

- $m$ represents somatic mutation features
- $c$ represents copy-number alteration features
- $s$ represents mutational signature features
- $a$ is age
- $z$ is sex

The classifier then outputs posterior probabilities across cancer types:

$$
\hat{p} = \mathrm{softmax}(f_\theta(x)).
$$

This setup is conceptually simple, but it is powerful because each feature family captures a different aspect of tumor biology:

- **Somatic mutations** capture recurrent driver patterns such as KRAS, EGFR, BRAF, and others.
- **Copy-number alterations** capture broad genomic structure and amplification or deletion events.
- **Mutational signatures** capture historical processes such as tobacco exposure or DNA repair defects.
- **Age and sex** inject prior clinical context without overwhelming the genomic signal.

The model then uses confidence thresholds to decide when a prediction is strong enough to be clinically interpretable.

---

## What Makes the Paper Strong

### 1. It uses clinically realistic data

This is not a model trained on idealized whole-genome research data. It uses **targeted NGS panels**, which is exactly the type of data many hospitals already generate. That makes the paper much more relevant to translational oncology.

### 2. It moves beyond pure accuracy

The article does report standard classifier performance, but it does not stop there. The authors show that OncoNPC predictions:

- align with **germline polygenic risk enrichment** for the predicted cancer types
- define CUP subgroups with **different survival outcomes**
- expand opportunities for **genomically guided therapy**

That last point is especially important. The abstract reports a **2.2-fold increase** in patients with CUP who could have received genomically guided therapies when OncoNPC predictions were taken into account.

### 3. It is interpretable enough to be clinically discussable

The paper uses SHAP-based feature interpretation rather than treating the model as a total black box. In the article's extended-data discussion, the authors show examples where mutational signatures, KRAS mutation, and copy-number features contribute to a specific tissue-of-origin prediction.

For clinical machine learning, this matters. A prediction that says "likely lung" is more convincing if it is accompanied by smoking-associated mutational signatures, lung-relevant driver patterns, and other coherent evidence.

{% include figure.liquid loading="eager" path="assets/img/blog/onconpc/figure4-shap-importance.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 4. SHAP-based feature importance for OncoNPC predictions. Somatic mutations (e.g., BRAF V600E, KRAS), mutational signatures (e.g., smoking, UV), copy-number alterations, and clinical variables each contribute to tissue-of-origin inference. Feature importance values reproduced from Moon et al., Nature Medicine 2023.
</div>

---

## Results That Matter Most

The abstract already contains the key numbers:

- trained on **36,445 tumors**
- **22** cancer types
- **weighted F1 = 0.942** for high-confidence held-out predictions
- high-confidence predictions in **41.2%** of **971** CUP tumors
- treatment-concordant cases showed improved outcomes with **HR = 0.348**
- **2.2-fold increase** in potentially genomically guided therapies

Those are unusually strong translational signals for a paper in this space.

{% include figure.liquid loading="eager" path="assets/img/blog/onconpc/figure3-confusion-matrix.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 3. Classification performance of OncoNPC across 22 cancer types (held-out set). The confusion matrix shows high accuracy along the diagonal for most cancer types, with the model achieving a weighted F1 score of 0.942 for high-confidence predictions. Data reproduced from Moon et al., Nature Medicine 2023.
</div>

Why? Because many papers on cancer classifiers stop at "the confusion matrix looks good." OncoNPC goes further and asks whether the classifier changes what you can _do_ for the patient. That is a much more demanding question.

The paper does not claim that every CUP case becomes solvable. It shows something more realistic and more useful: for a meaningful subset of cases, targeted sequencing can provide enough evidence to improve classification confidence and potentially guide therapy.

---

## Why This Is a Precision Medicine Paper, Not Just an AI Paper

Plenty of machine-learning papers in cancer are ultimately about prediction alone. OncoNPC is different because it fits the actual precision-medicine pipeline:

1. **Molecular profiling** of a patient tumor
2. **Probabilistic disease assignment**
3. **Therapy matching**
4. **Outcome association**

That structure is what makes the paper worth reading even if you are not interested in CUP specifically.

If your broader interest is precision medicine, this paper is a useful template for how to think:

- use data that already exist in clinical workflows
- predict something clinicians genuinely need
- connect the prediction to treatment options
- evaluate whether the prediction is associated with outcomes

That mindset is transferable to many other domains, including immunotherapy response prediction, molecular subtype assignment, and biomarker-guided patient stratification.

{% include figure.liquid loading="eager" path="assets/img/blog/onconpc/figure5-survival-concordance.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 5. Kaplan-Meier survival analysis of CUP patients stratified by treatment concordance with OncoNPC predictions. Patients receiving first palliative treatment concordant with the model's prediction showed significantly improved overall survival (HR = 0.348). Data reproduced from Moon et al., Nature Medicine 2023.
</div>

{% include figure.liquid loading="eager" path="assets/img/blog/onconpc/figure6-actionability.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 6. (A) OncoNPC expands genomically guided therapy opportunities by 2.2-fold for CUP patients. (B) Distribution of prediction confidence across 971 CUP tumors, with 41.2% receiving high-confidence classifications. Data reproduced from Moon et al., Nature Medicine 2023.
</div>

{% include figure.liquid loading="eager" path="assets/img/blog/onconpc/figure2-clinical-impact.svg" class="img-fluid rounded z-depth-1" zoomable=true %}

<div class="caption">
    Figure 2. Clinical interpretation of OncoNPC. The model is useful because it turns diagnostic uncertainty into a ranked, evidence-linked hypothesis that can influence treatment concordance and genomic actionability.
</div>

---

## A Useful Conceptual Model

One of the easiest ways to misunderstand OncoNPC is to think of it as a "truth machine" for tissue of origin. It is not. CUP is difficult precisely because ground truth is often unavailable.

A better way to think about the model is:

- it defines a **posterior over likely primary sites**
- it provides **confidence-aware prioritization**
- it offers a **decision-support layer** on top of pathology and sequencing

That framing is more realistic for clinical AI.

If a CUP tumor receives a high-confidence OncoNPC label of non-small cell lung cancer, supported by a smoking-associated mutational signature and lung-relevant genomic features, the classifier does not magically prove the origin. But it can raise the probability enough to influence treatment selection, molecular interpretation, and case review.

In other words, the model is valuable because it helps clinicians make **better bets under uncertainty**.

---

## Limitations and Caveats

This is an important paper, but it is not the last word.

### 1. The analysis is retrospective

The strongest clinical claim in the paper is based on retrospective observational analysis. That is informative, but it is not the same as a prospective trial where treatment is explicitly assigned using the model output.

### 2. The label space is constrained

OncoNPC models 22 cancer types. Real-world oncology is messier than that. Rare cancers, unusual histologies, and mixed phenotypes remain challenging.

### 3. Targeted panels are practical, but incomplete

Panel sequencing captures many clinically useful events, but it does not provide the full molecular context that whole-genome, methylation, transcriptomic, or histopathology-based approaches could add.

### 4. Concordance is not the same as causality

Patients receiving treatment concordant with model predictions did better, but retrospective concordance analyses can still be affected by selection bias, clinician judgment, performance status, and access to care.

### 5. Tissue-of-origin is only one axis

Precision oncology increasingly depends on more than origin alone. Molecular subtype, microenvironment, immune context, and lineage plasticity all matter. A future system will probably combine:

- genomics
- histopathology
- methylation
- RNA expression
- clinical notes

OncoNPC is strong precisely because it is practical, but the long-term direction is likely multimodal.

---

## What I Would Pair It With

If I were extending this paper in a modern research program, I would connect it to three additional layers:

### Single-cell and cellular-state context

Can the predicted tissue of origin also be supported by tumor microenvironment signatures or lineage-state markers?

### Regulatory interpretation

Can specific variants or copy-number changes identified by the model be linked to transcriptional consequences using expression or chromatin data?

### Actionability ranking

Instead of stopping at tissue prediction, can we move directly to:

- likely primary site
- likely pathway dependencies
- likely druggable targets
- likely immunotherapy relevance

That would make the pipeline even more aligned with real precision medicine.

---

## Key Takeaways

- OncoNPC addresses a real clinical bottleneck: assigning likely tissue of origin in cancers of unknown primary.
- The model uses targeted tumor sequencing plus basic clinical variables, which makes it practical for real oncology workflows.
- Its value is not just classifier accuracy. The important result is the connection to treatment concordance, survival, and actionable therapy matching.
- This paper is a strong example of what clinical machine learning should look like: existing assay, interpretable features, clinically meaningful output, and downstream decision relevance.
- For anyone interested in precision medicine, OncoNPC is less about CUP specifically and more about a translational pattern worth reusing.

---

## References and Resources

- Nature Medicine article: https://www.nature.com/articles/s41591-023-02482-6
- OncoNPC code repository: https://github.com/itmoon7/onconpc
