# LAMP RFLP

[![npm version](https://badge.fury.io/js/lamprflp.svg)](https://www.npmjs.com/package/lamprflp)

This is the JS module used by the site: http://creisle.github.io/creisle.lamprflp/

This computes the lamp products from a FASTA input file and outputs the fragments required
for the visual

## Getting Started

Install from npm

```bash
npm install lamprflp
```

then this can be used in a node module as follows

```js
const {createLampReaction} = require('lamprflp');

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
