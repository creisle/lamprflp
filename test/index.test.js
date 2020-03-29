const {
    parseLampFasta,
    createLampReaction,
    createAmplificationBlocks,
    getPrimerLayout,
    complementSeq
} = require('../src');

const data = require('./data.json');

const defaultFasta = `
> ${data.input.targetName}
${data.input.target}
> b1
${data.input.b1}
> b2
${data.input.b2}
> f1
${data.input.f1}
> f2
${data.input.f2}
`;

describe('parseLampFasta', () => {
    test('parses example ok', () => {
        const parsed = parseLampFasta(defaultFasta);
        expect(parsed).toEqual(data.input);
    });

    test('error on missing primer', () => {
        const fasta = `
> ${data.input.targetName}
${data.input.target}
> b1
${data.input.b1}
> b2
${data.input.b2}
> f1
${data.input.f1}
        `;
        expect(() => parseLampFasta(fasta)).toThrow('Missing sequence for f2 primer');
    });

    test('error on missing target', () => {
        const fasta = `
> b1
${data.input.b1}
> b2
${data.input.b2}
> f1
${data.input.f1}
> f2
${data.input.f2}
        `;
        expect(() => parseLampFasta(fasta)).toThrow('only primers detected, missing target sequence');
    });
});

describe('complementSeq', () => {
    test('complement only', () => {
        expect(complementSeq('atgc', false)).toEqual('tacg');
    });

    test('reverse and complement', () => {
        expect(complementSeq('atgc', true)).toEqual('gcat');
        expect(complementSeq('tgttgtgtggaattgtgagcggat', true)).toEqual('atccgctcacaattccacacaaca');
        expect(complementSeq('gtgcgggcctcttcgctattac', true)).toEqual('gtaatagcgaagaggcccgcac');
        expect(complementSeq('cgactctagaggatccccgggtac', true)).toEqual('gtacccggggatcctctagagtcg');
        expect(complementSeq('acaacgtcgtgactgggaaaaccct', true)).toEqual('agggttttcccagtcacgacgttgt');
    });
});

describe('getPrimerLayout', () => {
    test('matches known example', () => {
        const [reaction, layout] = getPrimerLayout(parseLampFasta(defaultFasta));
        expect(reaction).toEqual(data.reaction);
        expect(layout).toEqual(data.primerPositions);
    });
});


describe('createAmplificationBlocks', () => {
    test('matches known example', () => {
        const [reaction, layout] = getPrimerLayout(parseLampFasta(defaultFasta));
        const blocks = createAmplificationBlocks(reaction, layout);
        expect(blocks).toEqual(data.blocks);
    });
});

describe('createLampReaction', () => {
    let result = {};

    const {fragments} = data;

    beforeAll(() => {
        const [
            ladder100,
            ladder50,
            undigested,
            RsaI,
            BslI,
            HpaII,
            HinfI
        ] = createLampReaction(defaultFasta, ['100bp', '50bp', '', 'RsaI', 'BslI', 'HpaII', 'HinfI']);
        result = {ladder100, ladder50, undigested, RsaI, BslI, HpaII, HinfI};
    });

    test('undigested product ok', () => {
        expect(result.undigested).toEqual(fragments.undigested);
    });

    test.todo('100bp ladder ok');

    test.todo('50 bp ladder ok');

    test('RsaI product ok', () => {
        expect(result.RsaI).toEqual(fragments.RsaI);
    });

    test('BslI ok', () => {
        expect(result.BslI).toEqual(fragments.BslI);
    });

    test('HpaII ok', () => {
        expect(result.HpaII).toEqual(fragments.HpaII);
    });

    test('HinfI ok', () => {
        expect(result.HinfI).toEqual(fragments.HinfI);
    });

});