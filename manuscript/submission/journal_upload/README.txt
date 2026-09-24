Files for the JCAMD submission portal's "Figures and tables", "Supplementary
material", and "Related files" upload steps.

figures/
  Figure_1.tiff -- the single figure in the manuscript (native vs. redocked
  pose overlay for 9L40 and 9L4B), converted from the source PNG to TIFF.
  Source resolution is 1840x750 px; if the journal requests higher-
  resolution originals at the accepted-paper stage, this would need to be
  regenerated from the PyMOL session used to make the original figure
  (see the pipeline repository, validation_analysis/make_overlay_figure.pml),
  not upsampled from this file.

tables/
  Table_1.csv through Table_4.csv -- the four in-manuscript tables, each
  labeled and formatted to match the numbering and content in the
  manuscript body exactly.

supplementary/
  Supplementary_Information.pdf (and its .tex/.bbl source) -- combined
  supplementary document covering every "Supplementary Table/Section SN"
  referenced in the manuscript text (S1, S2, S4, S5, S6; there is no S3
  in the current manuscript). Two things to note before submitting:
    - Table S2 (the ~400-residue deletion list) could not be produced:
      no such enumerated list exists anywhere in the pipeline repository
      as a saved file, only the programmatic deletion step and its count.
      This is stated plainly in the document rather than filled in with
      invented residue IDs.
    - Table S4's note flags a small discrepancy: the main manuscript text
      says 680 of 786 derivatives passed the developability rules; the
      pipeline's current data file (developability_predictions.csv) shows
      677. This was not resolved -- see the note in the supplementary
      document and consider checking it before final submission.
  Table_S4/S5/S6 .csv -- the full underlying data for those three tables,
  provided alongside the typeset PDF versions (S4 in particular is long:
  677 rows, only the first 20 are typeset in the PDF for readability).

related_files/
  Cover_Letter.txt -- addressed to the JCAMD editor, matching the current
  manuscript's title, author list, corresponding author (Aryan Padarthi),
  and results exactly (the previous draft in the pipeline repository was
  written against an earlier title and an earlier corresponding author
  and was not reused as-is).
  Manuscript_Development_History.txt -- a plain-language summary of the
  five external review rounds this manuscript went through, plus the
  fabricated-content issue caught and corrected during final preparation
  (see the manuscript's own AI Assistance declaration). Offered for
  reviewer/editor context; not intended for publication.

What was not created: the portal's "Related files" description also
mentions "papers with overlapping authorship under consideration" -- none
of the authors have another such paper, so nothing was created for that.
"Lab validation reports" and "personal correspondence" beyond the review
history above did not apply either.
