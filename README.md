# LAMP RFLP

This is the JS module used by the site: http://creisle.github.io/creisle.lamprflp/

This computes the lamp products from a FASTA input file and outputs the fragments required
for the visual

## Getting Started

This module can be used as a node module as follows

```js
const {createLampReaction} = require('./src');

const {lanes} = createLampReaction(`
> b1
cgactctagaggatccccgggtac
> b2
tgttgtgtggaattgtgagcggat
> f1
acaacgtcgtgactgggaaaaccct
> f2
gtgcgggcctcttcgctattac
> seq
AATGCTACTACTATTAGTAG...TATGATTTATTGGATGTT
`, ['100bp', '50bp', '', 'RsaI']);
```

## Running the Tests

Tests are written with Jest and can be run as follows

```bash
npm test
```