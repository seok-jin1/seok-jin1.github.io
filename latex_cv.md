\documentclass[10pt, letterpaper]{article}

\usepackage{kotex}

% Packages:
\usepackage[
    ignoreheadfoot, % set margins without considering header and footer
    top=2 cm, % seperation between body and page edge from the top
    bottom=2 cm, % seperation between body and page edge from the bottom
    left=2 cm, % seperation between body and page edge from the left
    right=2 cm, % seperation between body and page edge from the right
    footskip=1.0 cm, % seperation between body and footer
    % showframe % for debugging 
]{geometry} % for adjusting page geometry
\usepackage{titlesec} % for customizing section titles
\usepackage{tabularx} % for making tables with fixed width columns
\usepackage{array} % tabularx requires this
\usepackage[dvipsnames]{xcolor} % for coloring text
\definecolor{primaryColor}{RGB}{0, 0, 0} % define primary color
% \usepackage{enumitem} % for customizing lists
\usepackage[revnum]{enumitem}
\usepackage{fontawesome5} % for using icons
\usepackage{amsmath} % for math
\usepackage[
    pdftitle={Seok-Jin Kang's CV},
    pdfauthor={Seok-Jin Kang},
    pdfcreator={LaTeX with RenderCV},
    colorlinks=true,
    urlcolor=cyan!60!blue,
    linkcolor=black,
]{hyperref} % for links, metadata and bookmarks
\usepackage[pscoord]{eso-pic} % for floating text on the page
\usepackage{calc} % for calculating lengths
\usepackage{bookmark} % for bookmarks
\usepackage{lastpage} % for getting the total number of pages
\usepackage{changepage} % for one column entries (adjustwidth environment)
\usepackage{paracol} % for two and three column entries
\usepackage{ifthen} % for conditional statements
\usepackage{needspace} % for avoiding page brake right after the section title
\usepackage{iftex} % check if engine is pdflatex, xetex or luatex

% Ensure that generate pdf is machine readable/ATS parsable:
\ifPDFTeX
    \input{glyphtounicode}
    \pdfgentounicode=1
    \usepackage[T1]{fontenc}
    \usepackage[utf8]{inputenc}
    \usepackage{lmodern}
\fi

\usepackage{charter}

% Some settings:
\raggedright
\AtBeginEnvironment{adjustwidth}{\partopsep0pt} % remove space before adjustwidth environment
\pagestyle{empty} % no header or footer
\setcounter{secnumdepth}{0} % no section numbering
\setlength{\parindent}{0pt} % no indentation
\setlength{\topskip}{0pt} % no top skip
\setlength{\columnsep}{0.15cm} % set column seperation
\pagenumbering{gobble} % no page numbering

%\titleformat{\section}{\needspace{4\baselineskip}\bfseries\large}{}{0pt}{}[\vspace{1pt}\titlerule]
\definecolor{sectionColor}{RGB}{0, 153, 153} 
\titleformat{\section}
  {\needspace{4\baselineskip}\color{sectionColor}\bfseries\large}{}{0pt}{}[\vspace{-6pt}\textcolor{black}{\titlerule}]

\titlespacing{\section}{
    % left space:
    -1pt
}{
    % top space:
    0.3 cm
}{
    % bottom space:
    0.2 cm
} % section title spacing

\renewcommand\labelitemi{$\vcenter{\hbox{\small$\bullet$}}$} % custom bullet points
\newenvironment{highlights}{
    \begin{itemize}[
        topsep=0.10 cm,
        parsep=0.10 cm,
        partopsep=0pt,
        itemsep=0pt,
        leftmargin=0 cm + 10pt
    ]
}{
    \end{itemize}
} % new environment for highlights


\newenvironment{highlightsforbulletentries}{
    \begin{itemize}[
        topsep=0.10 cm,
        parsep=0.10 cm,
        partopsep=0pt,
        itemsep=0pt,
        leftmargin=10pt
    ]
}{
    \end{itemize}
} % new environment for highlights for bullet entries

\newenvironment{onecolentry}{
    \begin{adjustwidth}{
        0 cm + 0.00001 cm
    }{
        0 cm + 0.00001 cm
    }
}{
    \end{adjustwidth}
} % new environment for one column entries

\newenvironment{twocolentry}[2][]{
    \onecolentry
    \def\secondColumn{#2}
    \setcolumnwidth{\fill, 4.5 cm}
    \begin{paracol}{2}
}{
    \switchcolumn \raggedleft \secondColumn
    \end{paracol}
    \endonecolentry
} % new environment for two column entries

\newenvironment{threecolentry}[3][]{
    \onecolentry
    \def\thirdColumn{#3}
    \setcolumnwidth{, \fill, 4.5 cm}
    \begin{paracol}{3}
    {\raggedright #2} \switchcolumn
}{
    \switchcolumn \raggedleft \thirdColumn
    \end{paracol}
    \endonecolentry
} % new environment for three column entries

\newenvironment{header}{
    \setlength{\topsep}{0pt}\par\kern\topsep\centering\linespread{1.5}
}{
    \par\kern\topsep
} % new environment for the header

\newcommand{\placelastupdatedtext}{% \placetextbox{<horizontal pos>}{<vertical pos>}{<stuff>}
  \AddToShipoutPictureFG*{% Add <stuff> to current page foreground
    \put(
        \LenToUnit{\paperwidth-2 cm-0 cm+0.05cm},
        \LenToUnit{\paperheight-1.0 cm}
    ){\vtop{{\null}\makebox[0pt][c]{
        \small\color{gray}\textit{Last updated in September 2024}\hspace{\widthof{Last updated in September 2024}}
    }}}%
  }%
}%

% save the original href command in a new command:
\let\hrefWithoutArrow\href

% new command for external links:


\begin{document}
    \newcommand{\AND}{\unskip
        \cleaders\copy\ANDbox\hskip\wd\ANDbox
        \ignorespaces
    }
    \newsavebox\ANDbox
    \sbox\ANDbox{$|$}

    \begin{header}
        \fontsize{25 pt}{25 pt}\selectfont Seok-Jin Kang

        \vspace{5 pt}

        \normalsize
        %\mbox{Biotechnology, Korea University}%
        %\kern 5.0 pt%
        %\AND%
        %\kern 5.0 pt%
        \mbox{\hrefWithoutArrow{mailto:laughingkang@korea.ac.kr}{laughingkang@korea.ac.kr}}%
        \kern 5.0 pt%
        \AND%
        \kern 5.0 pt%
        % \mbox{\hrefWithoutArrow{tel:+82-10-8258-1318}{+82-10-8258-1318}}%
        % \kern 5.0 pt%
        % \AND%
        % \kern 5.0 pt%
        \mbox{\textcolor{cyan!50!blue}{\hrefWithoutArrow{https://seok-jin1.github.io}{Github}}}%
        \kern 5.0 pt%
        \AND%
        \kern 5.0 pt%
        \mbox{\hrefWithoutArrow{https://scholar.google.com/citations?hl=ko&user=iexY3PEAAAAJ&view_op=list_works&sortby=pubdate}{Google Scholar}}%
        %\kern 5.0 pt%
        %\AND%
        %\kern 5.0 pt%
        %\mbox{\hrefWithoutArrow{https://linkedin.com/in/yourusername}{linkedin.com/in/yourusername}}%
        %\kern 5.0 pt%
        %\AND%
        %\kern 5.0 pt%
        %\mbox{\hrefWithoutArrow{https://github.com/yourusername}{github.com/yourusername}}%
    \end{header}

    \vspace{5 pt}

    \section{Bibliography}

    \begin{onecolentry}
        \textbf{Seok-Jin Kang, Ph.D.} is a self-motivated immunologist and computational biologist whose academic and research trajectory bridges traditional immunology and modern data-driven biology. During his M.S. and Ph.D. training in the \textit{Immune Modulation Laboratory (PI: Taehoon Chun)} at Korea University, he focused on understanding immune regulation, macrophage polarization, and the molecular mechanisms of host–pathogen interactions using both \textit{in vitro} and \textit{in vivo} models.
        
        \vspace{0.2cm}
        
        Recognizing that future scientific breakthroughs would emerge from integrating experimental and computational disciplines, Dr. Kang made a deliberate pivot from wet-lab immunology to data-driven research. Driven by curiosity and a strong ability to learn new concepts quickly, he has expanded his expertise beyond experimental immunology to computational biology, integrating \textbf{machine learning}, \textbf{deep learning}, and \textbf{quantum-inspired approaches} with \textbf{multi-omics} datasets. His recent publications explore the use of \textbf{AI-driven frameworks} for protein design, \textbf{TCR–peptide–MHC} interaction prediction, and the modeling of \textbf{intrinsically disordered regions (IDRs)} using \textbf{quantum neural networks}.
        
        \vspace{0.2cm}
        
        His long-term research goal is to advance the field of \textbf{precision medicine} by developing \textbf{multi-omics-based computational models} that can predict immune responses, guide therapeutic design, and enable personalized intervention strategies. As a fast learner who thrives in interdisciplinary environments, Dr. Kang seeks to bridge the gap between experimental immunology and computational systems biology to drive innovation in next-generation immunotherapies and biomedical engineering.
    \end{onecolentry}

    
    \section{Education}

    \begin{twocolentry}{
        \textbf{Mar 2018 – Feb 2025}
    }
        \textbf{Korea University} \\
        \textit{M.S. \& Ph.D. in Biotechnology } \\
        \textit{Immune Modulation Laboratory (Advisor: Prof. Taehoon Chun)}
    \end{twocolentry}
    
    \vspace{0.20 cm}
    
    \begin{twocolentry}{
        \textbf{Mar 2014 – Feb 2018}
    }
        \textbf{Korea University} \\
        \textit{B.S. in Biotechnology }
        
    \end{twocolentry}
    
    \vspace{0.10 cm}
    
    \section{Experience}
        
        \begin{twocolentry}{
    Mar 2025 – Present
}
        \textbf{Postdoctoral Fellow}, Korea University \\
        \textit{Immune Modulation Laboratory (PI: Taehoon Chun)}
        \end{twocolentry}
        \vspace{0.10 cm}
        \begin{onecolentry}
            \begin{highlights}
                \item Leading computational and experimental projects integrating \textbf{multi-omics} data with \textbf{machine learning and quantum-inspired models} to predict immunogenicity and design immune-specific biosensors.
                \item Applying \textbf{protein design frameworks} (AlphaFold, ProteinMPNN, and Rosetta) to engineer TCR–peptide–MHC interfaces for enhanced CAR-T specificity.
                \item Developing AI-driven pipelines for \textbf{CRISPR-based gene editing prediction} and analyzing \textbf{iPSC-derived T cell differentiation} using single-cell transcriptomics.
            \end{highlights}
        \end{onecolentry}

        \vspace{0.2 cm}
        



    \section{Publications}
    
    \begin{onecolentry}
        \textit{* indicates first author; † indicates principal investigator.}
    \end{onecolentry}
    
    \vspace{0.2cm}
    
    \begin{onecolentry}
        \textbf{Featured Publications}
        \begin{highlightsforbulletentries}
            \item \mbox{\textbf{\textit{SJ Kang}}*}, H Shin†, Biophysical mechanisms of spider-silk constituting element–induced stick-slip behavior and hydrogen bond regeneration for high toughness in silk fibers, \textit{International Journal of Biological Macromolecules}, 147027, 2025
    
            \item \mbox{\textbf{\textit{SJ Kang}}*}, H Shin†, Amino acid sequence-based IDR classification using ensemble machine learning and quantum neural networks, \textit{Computational Biology and Chemistry}, 108480, 2025
    
            \item \mbox{\textbf{\textit{SJ Kang}}*}, J Yang, NY Lee, CH Lee, IB Park, SW Park, HJ Lee, HW Park, HS Yun, T Chun†, Monitoring cellular immune responses after consumption of selected probiotics in immunocompromised mice, \textit{Food Science of Animal Resources} 42 (5), 903, 2022
    
            \item \mbox{\textbf{\textit{SJ Kang}}*}, IB Park, T Chun†, Open reading frame 5 protein of porcine circovirus type 2 induces RNF128 (GRAIL) which inhibits mRNA transcription of interferon-$\beta$ in porcine epithelial cells, \textit{Research in Veterinary Science} 140, 79–82, 2021
    
            \item \mbox{\textbf{\textit{SJ Kang}}*}, T Chun†, Structural heterogeneity of the mammalian polycomb repressor complex in immune regulation, \textit{Experimental \& Molecular Medicine} 52 (7), 1004–1015, 2020
        \end{highlightsforbulletentries}
    \end{onecolentry}
    
    \vspace{0.3cm}
    
    \begin{onecolentry}
        \textbf{Under Review}
        \begin{highlightsforbulletentries}
            \item \mbox{\textbf{\textit{SJ Kang}}*}, H Shin†, Quantum-Enhanced Transfer Learning for IDR Binding Partner Prediction, \textit{IEEE Transactions on Computational Biology and Bioinformatics}, under review (Round 1), 2025
        \end{highlightsforbulletentries}
    \end{onecolentry}
    
    \vspace{0.3cm}
    
    \begin{onecolentry}
        \textbf{Published Articles}
        \begin{highlightsforbulletentries}
            \item JH Yoon*, GB Yeon*, H Lee, H An, J Oh, \textbf{\textit{SJ Kang}}, IB Park, SA Lim, SS Hwang, DS Kim, JH Kim, T Chun†, YX Fu, J Bae, Tumor-targeted delivery of CXCL10 by mesenchymal stromal cells potentiates adoptive T cell therapy to treat solid tumors, \textit{Biomedicine \& Pharmacotherapy} 192, 118579, 2025
    
            \item CH Lee*, HJ Lee, SW Park, J Shin, \textbf{\textit{SJ Kang}}, IB Park, HK Kim, T Chun†, Mutational analysis of pig tissue factor pathway inhibitor $\alpha$ to increase anti-coagulation activity in pig-to-human xenotransplantation, \textit{Biotechnology Letters} 46 (4), 521–530, 2024
    
            \item SW Park*, IB Park*, \textbf{\textit{SJ Kang}}, J Bae, T Chun†, Interaction between host cell proteins and open reading frames of porcine circovirus type 2, \textit{Journal of Animal Science and Technology} 65 (4), 698, 2023
    
            \item SP Choi*, SW Park, \textbf{\textit{SJ Kang}}, SK Lim, MS Kwon, HJ Choi, T Chun†, Monitoring mRNA expression patterns in macrophages in response to two different strains of probiotics, \textit{Food Science of Animal Resources} 43 (4), 703, 2023
    
            \item CY Choi*, CH Lee, J Yang, \textbf{\textit{SJ Kang}}, IB Park, SW Park, NY Lee, HB Hwang, HS Yun, T Chun†, Efficacies of potential probiotic candidates isolated from traditional fermented Korean foods in stimulating immunoglobulin A secretion, \textit{Food Science of Animal Resources} 43 (2), 346, 2023
    
            \item SP Choi*, YC Choi, J Yang, CY Choi, CH Lee, \textbf{\textit{SJ Kang}}, IB Park, T Chun†, Monitoring mRNA transcription of genes involved in early pregnancy from endometrium and peripheral blood mononuclear cells of pregnant pigs with different parity, \textit{Reproduction in Domestic Animals} 53 (6), 1594–1599, 2018
    
            \item CY Choi*, YC Choi, IB Park, CH Lee, \textbf{\textit{SJ Kang}}, T Chun†, The ORF5 protein of porcine circovirus type 2 enhances viral replication by dampening type I interferon expression in porcine epithelial cells, \textit{Veterinary Microbiology} 226, 50–58, 2018
        \end{highlightsforbulletentries}
    \end{onecolentry}
    


    \section{Patents}

    \begin{onecolentry}
        \textbf{Registration}
    \end{onecolentry}
    
    \vspace{0.10 cm}
    
    \begin{onecolentry}
        \begin{highlightsforbulletentries}
            \item T Chun, J Yang, CH Lee, \textbf{\textit{SJ Kang}}, IB Park, Chimeric antigen receptor specifically binding to CD38 and use thereof, \textit{Korean Intellectual Property Office}, 10-2021-0025895, October 2021
    
            \item T Chun, J Yang, CH Lee, \textbf{\textit{SJ Kang}}, IB Park, Chimeric antigen receptor specifically binding to VLA-4 and use thereof, \textit{Korean Intellectual Property Office}, 10-2021-0025905, March 2021
        \end{highlightsforbulletentries}
    \end{onecolentry}
    
    \vspace{0.2 cm}
    
    \begin{onecolentry}
        \textbf{Pending}
    \end{onecolentry}
    
    \vspace{0.10 cm}
    
    \begin{onecolentry}
        \begin{highlightsforbulletentries}
            \item \textbf{\textit{SJ Kang}}, H Shin, Method and system for classifying intrinsically disordered proteins based on amino acid sequences using quantum neural networks, \textit{Korean Intellectual Property Office}, 10-2025-0078871, June 2025
    
            \item T Chun, IB Park, \textbf{\textit{SJ Kang}}, HJ Lee, Regulation of BCL11A expression by LDB1 transcriptional regulator and uses thereof, \textit{Korean Intellectual Property Office}, 10-2025-0010845, January 2025
    
            \item T Chun, IB Park, \textbf{\textit{SJ Kang}}, HJ Lee, Regulation of PLZF expression by LDB1 transcriptional regulator and uses thereof, \textit{Korean Intellectual Property Office}, 10-2025-0010844, January 2025
    
            \item HS Jin, T Chun, J Yang, CH Lee, NY Lee, IB Park, \textbf{\textit{SJ Kang}}, SW Park, HJ Lee, Chimeric antigen receptor specifically binding to CD138 and use thereof, \textit{Korean Intellectual Property Office}, 10-2023-0058232, May 2023
    
            \item HS Jin, T Chun, J Yang, CH Lee, NY Lee, IB Park, \textbf{\textit{SJ Kang}}, SW Park, HJ Lee, Recombinant protein which recognizes CADM1 and the composition comprising the same for treating cancer, \textit{Korean Intellectual Property Office}, 10-2023-0061710, May 2023
        \end{highlightsforbulletentries}
    \end{onecolentry}

    \section{Certifications}

        \begin{onecolentry}
            \begin{highlightsforbulletentries}
                \item \textbf{Big Data Analysis Engineer}, \textit{Korea Data Agency}, 2025, BAE-010002862
                
                \item \textbf{ADsP (Advanced Data Analytics Semi-Professional)}, \textit{Korea Data Agency}, 2024, ADsP-040004911
            \end{highlightsforbulletentries}
        \end{onecolentry}
    
    \section*{Research Projects}
    
    \begin{enumerate}[leftmargin=*,label=\textbf{\arabic*.}]
        \item \textbf{Q2431261} \hfill \textit{02/2025 -- 12/2025} \\
        \textbf{Title:} Safety evaluation of viral vector-based gene therapy with respect to non-specific genome insertion \\
        \textbf{Source:} Ministry of Food and Drug Safety \\
        \textbf{Role:} Research Assistant (PI: Taehoon Chun) \\
        \textbf{Description:} Conducted both \textbf{wet-lab and computational analyses} to evaluate vector integration safety. Prepared and optimized \textbf{whole-genome sequencing (WGS) libraries}, followed by data acquisition using \textbf{next-generation sequencing (NGS)}. Performed \textit{multi-omics integration and genomic correlation analyses} to identify and quantify off-target insertion events, examining their association with multiple factors such as \textbf{dose}, \textbf{administration route}, and \textbf{vector type} for comprehensive safety profiling.
        
    
        \item \textbf{Q2515581} \hfill \textit{05/2025 -- 12/2025} \\
        \textbf{Title:} Optimization of HLA gene editing for the development of allogeneic stem cell therapy \\
        \textbf{Source:} Korea National Institute of Health \\
        \textbf{Role:} Research Assistant (PI: Dongho Geum) \\
        \textbf{Description:} Designed and established \textbf{HLA-edited induced pluripotent stem cells (iPSCs)} for generating \textbf{hypoimmunogenic cell lines}. Conducted differentiation experiments to derive functional immune-evading cells, and performed \textbf{in vivo validation using mouse models} to assess immune compatibility and engraftment efficiency.

    
        \item \textbf{R2222761} \hfill \textit{09/2022 -- 02/2025} \\
        \textbf{Title:} Regulation of immune cell differentiation by chromatin 3D structural changes \\
        \textbf{Source:} National Research Foundation of Korea \\
        \textbf{Role:} Research Assistant (PI: Taehoon Chun) \\
        \textbf{Description:} Investigated how \textbf{chromatin 3D structural reorganization mediated by the chromatin looper LDB1} regulates transcriptional programs across multiple immune cell types. Conducted \textbf{multi-omics integration (RNA-seq, ATAC-seq, ChIP-seq)} to map enhancer–promoter interactions and applied \textbf{deep learning–based structural bioinformatics} to predict LDB1-associated chromatin loops driving immune cell differentiation.

    
        \item \textbf{R2212861} \hfill \textit{04/2022 -- 12/2023} \\
        \textbf{Title:} Development of universal SARS-CoV-2 and sarbecovirus vaccine candidate using NDV vector platform \\
        \textbf{Source:} Korea Health Industry Development Institute \\
        \textbf{Role:} Research Assistant (PI: Kisoon Kim) \\
        \textbf{Description:} Modeled viral antigen epitopes using \textbf{computational structural biology (AlphaFold)} to enhance cross-neutralizing immune responses. Conducted \textbf{in vitro} and \textbf{in vivo} immunogenicity evaluations using \textbf{animal models}, and analyzed vaccine-induced responses through \textbf{ELISA} and \textbf{flow cytometry}.

    
        \item \textbf{R2131151} \hfill \textit{02/2022 -- 12/2025} \\
        \textbf{Title:} Development of cellular immune evaluation technology for viral vector-based genetic vaccines \\
        \textbf{Source:} Ministry of Food and Drug Safety \\
        \textbf{Role:} Research Assistant (PI: Taehoon Chun) \\
        \textbf{Description:} 내용 수정해야함 Designed AI-assisted pipelines for predicting \textbf{TCR–peptide–MHC interactions} and validated cellular responses experimentally using \textbf{ELISPOT and flow cytometry}.
    
        \item \textbf{R1923081} \hfill \textit{09/2019 -- 02/2022} \\
        \textbf{Title:} Mechanistic study of macrophage polarization regulated by chromatin remodeling factor Phc2 \\
        \textbf{Source:} National Research Foundation of Korea \\
        \textbf{Role:} Research Assistant (PI: Taehoon Chun) \\
        \textbf{Description:} 내용 수정해야함 Performed transcriptomic profiling of macrophages under distinct polarization states. Applied \textbf{machine learning–based differential expression and pathway analysis} to identify Phc2-dependent networks.
    
        \item \textbf{R1727591} \hfill \textit{01/2018 -- 12/2020} \\
        \textbf{Title:} Development of DNA markers to enhance resistance against chronic wasting diseases \\
        \textbf{Source:} Rural Development Administration \\
        \textbf{Role:} Research Assistant (PI: Taehoon Chun) \\
        \textbf{Description:} Conducted \textbf{sequence analysis and SNP genotyping} to identify host resistance markers. Developed computational pipelines using \textbf{R and Python} for comparative genomics.
    
        \item \textbf{R1715822} \hfill \textit{03/2018 -- 06/2022} \\
        \textbf{Title:} Development of novel immune cell therapy targeting multiple myeloma \\
        \textbf{Source:} National Research Foundation of Korea \\
        \textbf{Role:} Research Assistant (PI: Taehoon Chun) \\
        \textbf{Description:} Designed and optimized \textbf{chimeric antigen receptor (CAR)} constructs through \textbf{protein engineering and molecular cloning}. Selected target-binding proteins based on affinity and specificity, and measured \textbf{binding kinetics (K\textsubscript{d})} using \textbf{surface plasmon resonance (SPR)} for functional validation.
    
        \item \textbf{Q2029461} \hfill \textit{01/2021 -- 06/2023} \\
        \textbf{Title:} Production and efficacy evaluation of novel recombinant immunotherapeutics \\
        \textbf{Source:} Immunological Designing Lab Co., Ltd. \\
        \textbf{Role:} Research Assistant (PI: Taehoon Chun) \\
        \textbf{Description:} Produced recombinant therapeutic proteins via \textbf{molecular cloning and protein purification}. Evaluated binding affinity using \textbf{SPR and ELISA} assays.
    
        \item \textbf{Q2015891} \hfill \textit{06/2020 -- 02/2021} \\
        \textbf{Title:} Screening of immune-enhancing probiotics using immunodeficient mouse models \\
        \textbf{Source:} CJ CheilJedang \& BLOSSOM PARK \\
        \textbf{Role:} Research Assistant (PI: Taehoon Chun) \\
        \textbf{Description:} Conducted animal experiments and analyzed cytokine profiles via \textbf{flow cytometry and qPCR}. Applied \textbf{statistical analysis in R} to evaluate immune enhancement.
    
        \item \textbf{Q1815701} \hfill \textit{06/2018 -- 11/2018} \\
        \textbf{Title:} Determination of dosing period and amount of probiotics for alleviating allergic rhinitis \\
        \textbf{Source:} CJ CheilJedang \\
        \textbf{Role:} Research Assistant (PI: Taehoon Chun) \\
        \textbf{Description:} Evaluated probiotic-induced modulation of Th1/Th2 cytokine ratios using \textbf{ELISA and gene expression profiling}. Applied \textbf{statistical modeling} to identify optimal dosing.
    
        \item \textbf{PJ013271012018} \hfill \textit{03/2018 -- 12/2020} \\
        \textbf{Title:} Development of genetic markers associated to the resistance to PMWS \\
        \textbf{Source:} Rural Development Administration \\
        \textbf{Role:} Research Assistant (PI: Taehoon Chun) \\
        \textbf{Description:} Performed comparative genomics and \textbf{sequence-based variant analysis}. Utilized \textbf{bioinformatics tools (Biopython, BLAST)} to detect resistance-associated mutations.
    
        \item \textbf{Q1613841} \hfill \textit{07/2016 -- 04/2017} \\
        \textbf{Title:} Selection and in vivo efficacy verification of probiotics with enhanced IgA secretion \\
        \textbf{Source:} CJ CheilJedang \\
        \textbf{Role:} Research Assistant (PI: Taehoon Chun) \\
        \textbf{Description:} Measured mucosal IgA levels in mouse models and identified strains with strong immunostimulatory potential. Applied \textbf{bioinformatics analysis of gene expression} and immune pathway mapping.
    \end{enumerate}
    

    
   


    \section{Book chapters}

    \begin{onecolentry}
        \textbf{International Publications}
        \begin{highlightsforbulletentries}
            \item \mbox{\textbf{\textit{SJ Kang}}}, \textbf{2026 Science Trends: Multi-omics}, \textit{Kindle}, ISBN 9798266588004, September 2025
    
            \item \mbox{\textbf{\textit{SJ Kang}}}, \textbf{2025 Science Trends: Artificial Intelligence}, \textit{Kindle}, ISBN 979-8307079270, January 2025
        \end{highlightsforbulletentries}
    
        \vspace{0.2cm}
    
        \textbf{Domestic Publications}
        \begin{highlightsforbulletentries}
            \item \mbox{\textbf{\textit{SJ Kang}}}, \textbf{2026 생명과학 트렌드: 멀티오믹스(2026 Life Science Trends: Multi-omics)}, \textit{Pubple}, ISBN 9788924169034, August 2025
    
            \item \mbox{\textbf{\textit{SJ Kang}}}, \textbf{2026 과학 트렌드: 개인 맞춤형 치료(2026 Science Trends: 개인 Personalized Medicine)}, \textit{Pubple}, ISBN 9788924169041, August 2025
    
            \item \mbox{\textbf{\textit{SJ Kang}}}, \textbf{2025 과학 트렌드: 인공지능(2025 Science Trends: Artificial Intelligence)}, \textit{Pubple}, ISBN 9788924143621, December 2024
        \end{highlightsforbulletentries}
    \end{onecolentry}

    \section{Research Interests}

    \begin{onecolentry}
        Precision Medicine, CAR-T Cell Therapy, Multi-omics Integration, Computer-aided Drug Design (CADD), Quantum Computing
    \end{onecolentry}
    

\end{document}