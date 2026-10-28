# PriorR 4.0 user guide

`guide.html` is the source of `app/www/PriorR4_UserGuide.pdf` (linked from the Help tab). The screenshots
(`img/`) were taken from PriorR with the synthetic demo data (`app/Program_files/examples/demo_*_synthetic.tsv`,
made by `app/Program_files/scripts/make_demo_data.R`), never with patient data.

To make the PDF again (Chrome or Chromium):

    google-chrome --headless --no-pdf-header-footer --print-to-pdf=app/www/PriorR4_UserGuide.pdf docs/user_guide/guide.html
