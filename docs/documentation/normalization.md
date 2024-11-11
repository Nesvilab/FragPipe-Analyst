# Normalization

By default, the input data are not normalized, as they are assumed to have already been normalized by the FragPipe quantification tools.
Optionally, normalization methods listed below could be applied:

- Variance-stabilizing normalization: performed using the R package vsn. It's available only for DDA- and DIA-based LFQ data.
- Median centered: quantification values of each samples are minused against median values of each sample.
