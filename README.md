# 3DNA (v2.4) — The Classic 3DNA Suite

[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](LICENSE)
[![Status](https://img.shields.io/badge/Status-Maintenance_Only-yellow.svg)]()
[![Simple Parameters](https://img.shields.io/badge/Simple_BP_Parameters-PDF-orange.svg)](CCN201601-base-pair-geometry.pdf)

**3DNA (v2.4)** is the final, comprehensive release of the classic **3DNA** suite of programs for the analysis, rebuilding, and visualization of three-dimensional nucleic acid structures.

> **Maintenance Status & Modern Successors:**
> **3DNA v2.4** is still maintained for bug fixes where appropriate, but **no new features are planned**. Existing functionalities of 3DNA have been superseded by [**DSSR**](http://x3dna.org/) and the [**wDSSR** web application](https://web.x3dna-dssr.org). In short, **3DNA represents the past, whereas DSSR and wDSSR represent the future.**

---

## Overview & Key Capabilities

Written in strict ANSI C with connecting Ruby (and legacy Perl) scripts, **3DNA v2.4** provides a self-contained command-line environment for structural bioinformatics of DNA, RNA, and their complexes:

* **Rigorous Local Base-Pair, Step, and Helical Parameters:** Uses the standard base-centered reference frame (Olson et al., 2001) and the matrix-based CEHS scheme (El Hassan & Calladine, 1995; Lu et al., 1997) to compute six local base-pair parameters (*Shear, Stretch, Stagger, Buckle, Propeller, Opening*), six dimeric step parameters (*Shift, Slide, Rise, Tilt, Roll, Twist*), and six local helical parameters (*x-displacement, y-displacement, Helical rise, Inclination, Tip, Helical twist*).
* **Intuitive "Simple" Base-Pair & Step Parameters:** In addition to mathematically rigorous local parameters (which are required for exact 3D rebuilding), 3DNA v2.4 computes six **simple base-pair parameters** and six **simple step parameters** (Lu & Olson, 2016; Li et al., 2019). Defined using torsion angles across the RC8–YC6 (or RN9–YN1) axis, simple parameters provide an intuitive qualitative characterization of non-Watson–Crick pairs (such as Hoogsteen pairs and G-tetrads) and their stacking steps. See **[`CCN201601-base-pair-geometry.pdf`](CCN201601-base-pair-geometry.pdf)** for details.
* **Exact 3D Structure Rebuilding & Fiber Models:** Reconstructs full-atom PDB structures (with or without sugar-phosphate backbones) and Calladine–Drew block schematics (`rebuild`), and generates 55 regular DNA and RNA fiber diffraction models (`fiber`).
* **Automated Base-Pair, Multiplet & Helix Detection (`find_pair`):** Identifies canonical and non-canonical base pairs, higher-order coplanar base multiplets (triplets, quadruplexes, pentaplets), and continuous or coaxial double-helical segments directly from PDB coordinates.
* **Ensemble Analysis & Utilities (`x3dna_ensemble`, `x3dna_utils`):** Ruby-based command-line workflows for analyzing NMR ensembles or MD trajectories and performing common coordinate manipulations (such as base frame orientations and *in silico* base mutations).
* **Reproducible *Nature Protocols* Recipes (`np_recipes/`):** All scripts, data files, and step-by-step protocols reported in Lu & Olson (2008) are included in the `np_recipes/` directory.

For a foundational tutorial and step-by-step mathematical walkthrough of the quaternion least-squares fitting, reference frame construction, and local step/helical rebuilding equations, please see the user manual in the **[`x3dna-v1.5` repository](https://github.com/x3dna-dssr/x3dna-v1.5)**.

---

## Directory Structure & Setup

### Package Organization

* **`config/`** — Standard base residue geometries (`Atomic_?.pdb` and `atomic/`), Calladine–Drew block definitions (`block/`), residue mapping (`baselist.dat`), and plotting parameter files (`misc_3dna.par`, `ps_image.par`, `fig_image.par`, `raster3d.par`).
* **`fiber/`** — Repeating unit coordinates (`str01`–`str55`) for 55 regular DNA and RNA fiber diffraction models.
* **`np_recipes/`** — Complete scripts (recipes) and associated data files to reproduce the 2008 *Nature Protocols* paper. *(Note: Generating raster graphics images requires external tools such as MolScript, Raster3D, PyMOL, ImageMagick, and Ghostscript; Xfig is useful for `.fig` files, and RasMol v2.6.x with `-alchemy -noconnect` for `.alc` files.)*
* **`examples/`** — Additional worked examples for analysis, rebuilding, NMR ensembles, Calladine–Drew schematics, stacking diagrams, and triplexes.
* **`src/`** — Complete ANSI C source code and `Makefile` for compiling the C executables into `bin/`.
* **`bin/` & `lib/`** — Core Ruby scripts (`blocview`, `x3dna_ensemble`, `x3dna_setup`, `x3dna_utils`) backed by `lib/*.rb` (`ensemble.rb`, `miscs.rb`, `parser.rb`, `trollop.rb`), plus the 15 compiled C executables once built (`alc2img`, `analyze`, `cehs`, `comb_str`, `ex_str`, `fiber`, `find_pair`, `frame_mol`, `get_part`, `pdb2img`, `r3d_atom`, `rebuild`, `rotate_mol`, `stack2img`, `std_base`).
* **`perl_scripts/`** — Original Perl scripts from earlier versions (v1.5 and v2.0) kept for backward reference (as of v2.1, 3DNA transitioned from Perl to Ruby, consolidating utilities into `x3dna_utils` and `x3dna_ensemble`).

### Compilation & Environment Configuration

1. **Compile the C executables:**
   ```bash
   cd src
   make
   cd ..
   ```

2. **Configure your shell environment:**
   Run `./bin/x3dna_setup` for interactive setup instructions, or add the `X3DNA` environment variable and `$X3DNA/bin` to your shell configuration file:

   **For `bash` / `zsh` (add to `~/.bashrc` or `~/.zshrc`):**
   ```bash
   export X3DNA=/path/to/x3dna-v2.4
   export PATH=$X3DNA/bin:$PATH
   ```

   **For `csh` / `tcsh` (add to `~/.cshrc`):**
   ```csh
   setenv X3DNA /path/to/x3dna-v2.4
   set path = ($path$X3DNA/bin)
   ```

---

## Citation

If you use or reference 3DNA in your work, please cite at least one of the primary 3DNA publications (1–2), along with the relevant methodological references (3–6) as appropriate:

### Primary 3DNA Publications

1. **3DNA *Nature Protocols* Paper:**
   Lu, X.-J., & Olson, W. K. (2008). 3DNA: a versatile, integrated software system for the analysis, rebuilding and visualization of three-dimensional nucleic-acid structures. *Nature Protocols*, 3(7), 1213–1227. [https://doi.org/10.1038/nprot.2008.104](https://doi.org/10.1038/nprot.2008.104)

2. **Original 3DNA *NAR* Paper:**
   Lu, X.-J., & Olson, W. K. (2003). 3DNA: a software package for the analysis, rebuilding and visualization of three-dimensional nucleic acid structures. *Nucleic Acids Research*, 31(17), 5108–5121. [https://doi.org/10.1093/nar/gkg680](https://doi.org/10.1093/nar/gkg680)

### Additional Methodological References

3. **Analyzing & Building Structures (*JoVE*):**
   Colasanti, A. V., Lu, X.-J., & Olson, W. K. (2013). Analyzing and building nucleic acid structures with 3DNA. *Journal of Visualized Experiments (JoVE)*, 74, e4401. [https://doi.org/10.3791/4401](https://doi.org/10.3791/4401)

4. **"Simple" Base-Pair Parameters (*CCN*):**
   Lu, X.-J., & Olson, W. K. (2016). Characterization of base pair geometry. *Computational Crystallography Newsletter*, 7(1), 6–9. ([PDF](CCN201601-base-pair-geometry.pdf))

5. **3DNA v2.4 Updates & Simple Step Parameters (*NAR*):**
   Li, S., Olson, W. K., & Lu, X.-J. (2019). Web 3DNA 2.0 for the analysis, visualization, and modeling of 3D nucleic acid structures. *Nucleic Acids Research*, 47(W1), W26–W34. [https://doi.org/10.1093/nar/gkz394](https://doi.org/10.1093/nar/gkz394)

6. **Major & Minor Groove Width Calculation:**
   El Hassan, M. A., & Calladine, C. R. (1998). Two distinct modes of protein-induced bending in DNA (Appendix: Calculation of Major- and Minor-groove Widths). *Journal of Molecular Biology*, 282(2), 331–343. [https://doi.org/10.1006/jmbi.1998.1994](https://doi.org/10.1006/jmbi.1998.1994)

### BibTeX

```bibtex
@article{Lu2008_3DNA_NatProtoc,
  author  = {Lu, Xiang-Jun and Olson, Wilma K.},
  title   = {{3DNA}: a versatile, integrated software system for the analysis, rebuilding and visualization of three-dimensional nucleic-acid structures},
  journal = {Nature Protocols},
  volume  = {3},
  number  = {7},
  pages   = {1213--1227},
  year    = {2008},
  doi     = {10.1038/nprot.2008.104}
}

@article{Lu2003_3DNA,
  author  = {Lu, Xiang-Jun and Olson, Wilma K.},
  title   = {{3DNA}: a software package for the analysis, rebuilding and visualization of three-dimensional nucleic acid structures},
  journal = {Nucleic Acids Research},
  volume  = {31},
  number  = {17},
  pages   = {5108--5121},
  year    = {2003},
  doi     = {10.1093/nar/gkg680}
}

@article{Colasanti2013_3DNA_JoVE,
  author  = {Colasanti, Andrew V. and Lu, Xiang-Jun and Olson, Wilma K.},
  title   = {Analyzing and building nucleic acid structures with {3DNA}},
  journal = {Journal of Visualized Experiments},
  number  = {74},
  pages   = {e4401},
  year    = {2013},
  doi     = {10.3791/4401}
}

@article{Lu2016_SimpleBP,
  author  = {Lu, Xiang-Jun and Olson, Wilma K.},
  title   = {Characterization of base pair geometry},
  journal = {Computational Crystallography Newsletter},
  volume  = {7},
  number  = {1},
  pages   = {6--9},
  year    = {2016}
}

@article{Li2019_w3DNA2,
  author  = {Li, Shuxiang and Olson, Wilma K. and Lu, Xiang-Jun},
  title   = {Web {3DNA} 2.0 for the analysis, visualization, and modeling of {3D} nucleic acid structures},
  journal = {Nucleic Acids Research},
  volume  = {47},
  number  = {W1},
  pages   = {W26--W34},
  year    = {2019},
  doi     = {10.1093/nar/gkz394}
}

@article{ElHassan1998_GrooveWidth,
  author  = {El Hassan, M. A. and Calladine, Christopher R.},
  title   = {Two distinct modes of protein-induced bending in {DNA}},
  journal = {Journal of Molecular Biology},
  volume  = {282},
  number  = {2},
  pages   = {331--343},
  year    = {1998},
  doi     = {10.1006/jmbi.1998.1994}
}
```

---

## Community & Support

For questions, comments, or discussions regarding the 3DNA/DSSR ecosystem, please visit [GitHub Discussions](https://github.com/orgs/x3dna-dssr/discussions).

## License

Distributed under the [BSD 3-Clause License](LICENSE).
