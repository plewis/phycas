from phycas import Newick
from phycas.Utilities.PDFTree import PDFTree

pdftree = PDFTree()
pdftree.pdf_newick                = Newick("('P. fimbriata':0.1,('P. articulata':0.09,'P. parksii':0.04)v:0.1,('P. macrophylla':0.14,(('P. gracilis':0.08,('P. ciliata':0.02,'P. basiramia':0.03)x:0.05)y:0.01,'P. polygama':0.09)z:0.06)w:0.07)u", Newick.TAXA_NAMES)
pdftree.pdf_outgroup_taxon        = 'P. macrophylla' # use None to display tree as is
pdftree.pdf_tip_label_font        = 'Times-Italic'   # must be one of the 14 standard font names
pdftree.pdf_tip_label_height      = 24               # in units of points
pdftree.pdf_scalebar_label_font   = 'Times-Roman'    # must be one of the 14 standard font names
pdftree.pdf_scalebar_label_height = 10               # in units of points
pdftree.pdf_ladderize             = 'right'          # valid choices: None, 'left', or 'right'
pdftree.pdf_scalebar_position     = 'bottom'         # value choices: None, 'top', or 'bottom'
pdftree.pdf_page_width            = 8.5              # in units of inches
pdftree.pdf_page_height           = 11.0             # in units of inches
pdftree.pdf_line_width            = 2.0              # in units of points
pdftree.pdf_left_margin           = 1.0              # in units of inches
pdftree.pdf_right_margin          = 1.0              # in units of inches
pdftree.pdf_top_margin            = 1.0              # in units of inches
pdftree.pdf_bottom_margin         = 1.0              # in units of inches
pdftree.pdf_filename              = 'test.pdf'
pdftree.pdftree()

