@default_files = ("report.tex");
$pdflatex = 'lualatex --file-line-error %O %S'; # use xelatex (faster) or lualatex (better supported)
$pdf_mode = 1;
