# GeneMap

GeneMap, A python tool for identifying the mapping of genes between two genomes by multi evidences.
![Flow chart](https://raw.githubusercontent.com/wjwei-handsome/wwjPic/main/img/20220722001234.png)

## Table of Contents

- [Background](#background)
- [Install](#install)
- [Usage](#usage)
- [Badge](#badge)
- [Maintainers](#maintainers)
- [Contributing](#contributing)
- [License](#license)

## Background

When we have multiple genomes of the same species (or related species), we usually need to identify the gene mapping relationship between the two genomes(Gene A in genome a, mapped to gene B in genome B). We used four alternative evidence, including:

1. RBH
   **R**eciprocal **B**est BLAST **H**its

   > Despite these and other complication, the identification of reciprocal best hits for gene products is a good first approximation to the identification of orthologues in two or more organisms. It forms the basis for many orthology-finding tools, such as MCL, OrthoMCL and OrthoFinder. It can be carried out by carrying out BLAST+ searches using a short program, and this will be illustrated below.

2. ortholog
   [OrthoFinder](https://github.com/davidemms/OrthoFinder)
3. synteny block
   [McScanX](<https://github.com/tanghaibao/jcvi/wiki/MCscan-(Python-version)>)
4. crosssmap
   [minimap](https://github.com/lh3/minimap2)
   or
   [AnchorWave](https://github.com/baoxingsong/AnchorWave)

## Install

TODO

```sh
$ pip3 install GeneMap
```

## Usage

```sh
$ GeneMap.py -l list.txt -d work_dir -q GenomeA -t GenomeB -o output_prefix
```

## Maintainers

[@wjwei-handsome](https://github.com/wjwei-handsome)

## Contributing

Feel free to dive in! [Open an issue](https://github.com/wjwei-handsome/GeneMap/issues/new) or submit PRs.

Standard Readme follows the [Contributor Covenant](http://contributor-covenant.org/version/1/3/0/) Code of Conduct.

### Contributors

Thanks [@songtaogui](https://github.com/songtaogui) for design ideas.

## License

[GPL-3.0](LICENSE) © Weiwenjie
