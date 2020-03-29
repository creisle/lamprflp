const {enzymes, bpLadders} = require('./data.json');

const PRIMERS = ['b1', 'b2', 'f1', 'f2', 'linker'];


/**
 * @typedef {object} LampReaction
 *
 * @property {string} b1 the b1 primer sequence
 * @property {string} b2 the b2 primer sequence
 * @property {string} f1 the f1 primer sequence
 * @property {string} f2 the f2 primer sequence
 * @property {string} linker the linker sequence
 * @property {string} target the target sequence
 * @property {string} targetName the name of the target sequence
 */

/**
 * @typedef {object} PrimerPositions
 *
 * @property {object} b1 the b1 primer position
 * @property {object} b2 the b2 primer position
 * @property {object} f1 the f1 primer position
 * @property {object} f2 the f2 primer position
 */

/**
 * @typedef {object} LampBuildingBlocks
 *
 * @property {string} bplus
 * @property {string} bminus
 * @property {string} fplus
 * @property {string} fminus
 * @property {string} plus
 * @property {string} minus
 */

/**
 *
 * @param {string} fasta  the input fasta content
 *
 * @returns {LampReaction} the parsed reaction parameters
 */
const parseLampFasta = (fasta) => {

    const sequences = {linker: 'tttt'};

    const lines = fasta.split(/\n/g);
    let state = '';

    for (const rawLine of lines) {
        const line = rawLine.trim();

        if (!line) {
            continue;
        }

        if (!state) {
            if (!line.startsWith('>')) {
                throw new Error(`expected a new sequence name marker '>' but found ${line}`);
            }
            state = line.replace('>', '').trim();
            sequences[state] = '';
        } else {
            if (line.startsWith('>')) {
                if (!sequences[state]) {
                    throw new Error(`expected sequence content for marker ${state} but found ${line}`);
                } else {
                    state = line.replace('>', '').trim();
                    sequences[state] = '';
                }
            } else {
                sequences[state] += line;
            }
        }
    }
    const reaction = {};

    for (const primer of PRIMERS) {
        if (!sequences[primer]) {
            throw new Error(`Missing sequence for ${primer} primer`);
        }
        reaction[primer] = sequences[primer].toLowerCase();
    }

    const targets = Object.keys(sequences).filter(seq => !PRIMERS.includes(seq));

    if (!targets.length) {
        throw new Error('only primers detected, missing target sequence');
    } else if (targets.length !== 1) {
        throw new Error(`too many targets detected (${targets.join(', ')}). Only 1 target sequence is supported`);
    }


    [reaction.targetName] = targets;
    reaction.target = sequences[reaction.targetName].toLowerCase();
    return reaction;
};

const complementSeq = (seq, reverse = false) => {
    if (seq == null) {
        return null;
    }
    seq = seq.toLowerCase();
    var result = [];

    for (var k = 0; k < seq.length; k++){
        var i = k;

        if (reverse){
            i = seq.length - k - 1;
        }

        switch (seq.charAt(i)){
            case 'c':
                result.push('g');
                break;
            case 'g':
                result.push('c');
                break;
            case 't':
                result.push('a');
                break;
            case 'a':
                result.push('t');
                break;
            case 'r': //a, g --> t, c
                result.push('y');
                break;
            case 'y': //c, t --> g, a
                result.push('r');
                break;
            case 'k': //g, t --> c, a
                result.push('m');
                break;
            case 'm': //a, c --> t, g
                result.push('k');
                break;
            case 'b': //c, g, t --> g, c, a
                result.push('v');
                break;
            case 'd': //a, g, t --> t, c, a
                result.push('h');
                break;
            case 'h': //a, c, t --> t, g, a
                result.push('d');
                break;
            case 'v': //a, c, g --> t, g, c
                result.push('b');
                break;
            default: //s, w, n
                result.push(seq.charAt(i));
                break;
        }
    }
    return result.join('');
};

/**
 * find the best alignment position for an primer on a target
 *
 * @param {string} targetSeq the target string
 * @param {string} primerSeq the primer string (primerSeq.length<str.length)
 *
 * @returns returns a match object with the index position and whether or not the best match was found with the original primer string or the reverse complement
 */
const alignPrimer = (targetSeq, primerSeq) => {
    var revcomp = complementSeq(primerSeq, true);
    var minscore = primerSeq.length; //min number of differences btwn primer and sequence
    var reverse = false; //min number of differences btwn revcomp and sequence
    var index = 0;

    for (var i = 0; i < targetSeq.length - primerSeq.length; i++){
        let score = 0,
            revscore = 0;

        for (var j = 0; j < primerSeq.length; j++){
            if (targetSeq.charAt(i + j) != primerSeq.charAt(j)){
                score += 1;
            }
            if (targetSeq.charAt(i + j) != revcomp.charAt(j)){
                revscore += 1;
            }
            if (score > minscore && revscore > minscore){
                break;
            }
        }

        if (score <= revscore && score < minscore){
            minscore = score; reverse = false; index = i;
        } else if (revscore < score && revscore < minscore){
            minscore = revscore; reverse = true; index = i;
        }
    }
    return {index: index, reverse: reverse};
};

/**
 * Check that the layout of the primers is valid and the reaction can proceed. Fix common mistakes
 * where it is possible to guess the correct input
 *
 * @param {LampReaction} reaction the parse lamp reaction paramters
 *
 * @returns {Array[LampReaction, PrimerPositions]} the modified lamp reaction and the primer matches
 */
const getPrimerLayout = (inputReaction) => {
    const rxn = Object.assign({}, inputReaction);
    var b2match = alignPrimer(rxn.target, rxn.b2);
    var f2match = alignPrimer(rxn.target, rxn.f2);

    if (b2match.index < f2match.index){
        if (b2match.reverse == true){
            console.warn('input string was the other strand from what was needed. will take the compliment and re-compute');
            rxn.target = complementSeq(rxn.target, false); //just gives the compliment, not reversed
            b2match = alignPrimer(rxn.target, rxn.b2);
            f2match = alignPrimer(rxn.target, rxn.f2);
        }
    } else {
        if (b2match.reverse == true){
            console.warn('need reverse compliment of the input strand. generating...');
            rxn.target = complementSeq(rxn.target, true); //just gives the compliment, not reversed
            b2match = alignPrimer(rxn.target, rxn.b2);
            f2match = alignPrimer(rxn.target, rxn.f2);
        } else {
            console.warn('just need to reverse the input strand');
            rxn.target = rxn.target.split('').reverse().join('');
            b2match = alignPrimer(rxn.target, rxn.b2);
            f2match = alignPrimer(rxn.target, rxn.f2);
        }
    }

    //after this all should be b2.....f2c

    var temp = rxn.target.substring(
        b2match.index + rxn.b2.length,
        f2match.index
    );

    if (temp.length < rxn.b1.length){
        throw new Error('error in the primer input. primers do not align to the target in the expected format. please check the provided sequences');
    }
    var b1match = alignPrimer(temp, rxn.b1);
    b1match.index += b2match.index + rxn.b2.length;
    temp = rxn.target.substring(b1match.index + rxn.b1.length, f2match.index);

    if (temp.length < rxn.b1.length){
        throw new Error('error in the primer input. primers do not align to the target in the expected format. please check the provided sequences');
    }
    var f1match = alignPrimer(temp, rxn.f1);
    f1match.index += b1match.index + rxn.b1.length;

    //all the indicies must be right due to our substring restriction. now check the directions
    //already know that b2 is not reversed and order is b2...b1....f1....f2
    if (!b1match.reverse){ //should be true, revcomp the primer
        rxn.b1 = complementSeq(rxn.b1, true);
        b1match.reverse = true;
    }
    if (f1match.reverse){ //should be false
        rxn.f1 = complementSeq(rxn.f1, true);
        f1match.reverse = true;
    }
    if (!f2match.reverse){ //should be true
        rxn.f2 = complementSeq(rxn.f2, true);
        f2match.reverse = true;
    }

    var matches = {b1: b1match, b2: b2match, f1: f1match, f2: f2match};
    return [rxn, matches];
};

/**
 * Create the building block sequences for the LAMP reaction
 * @param {LampReaction} reaction
 * @param {PrimerPosition} primerPositions
 *
 * @returns {LampBuildingBlocks} the building blocks for the amplification products
 */
const createAmplificationBlocks = (reaction, primerPositions) => {
    const bplus = `${reaction.b1}${reaction.linker}${reaction.b2}${
        reaction.target.substring(
            primerPositions.b2.index + reaction.b2.length,
            primerPositions.b1.index + reaction.b1.length
        )}`;
    const bminus = complementSeq(bplus, true);
    const fplus = `${
        reaction.target.substring(
            primerPositions.f1.index,
            primerPositions.f2.index - 1
        )}${reaction.f2}${reaction.linker}${reaction.f1}`;
    const fminus = complementSeq(fplus, true);
    const plus = reaction.target.substring(
        primerPositions.b1.index + reaction.b1.length,
        primerPositions.f1.index - 1
    );
    const minus = complementSeq(plus, true);
    return {bplus, bminus, fplus, fminus, plus, minus};
};

/**
 *
 * @param {LampBuildingBlocks} blocks
 */
const amplifyBlocks = ({bplus, bminus, plus, minus, fplus, fminus}) => {
    const products = [];
    let temp = `${bplus}${plus}${fplus}`; //6
    products.push(temp);
    temp += `${minus}${bminus}`; //7
    products.push(temp);
    temp += `${plus}${fminus}${minus}${bminus}`; //9
    products.push(temp);
    temp += `${plus}${fplus}${minus}${bplus}${plus}${fminus}${minus}${bminus}`; //13
    products.push(temp);
    temp += `${plus}${fplus}${minus}${bminus}${plus}${fminus}${minus}${bplus}${plus}${fplus}${minus}${bplus}${plus}${fminus}${minus}${bminus}`; //15
    products.push(temp);
    temp = `${fminus}${minus}${bplus}${plus}${fplus}${minus}${bplus}${plus}${fminus}${minus}${bminus}`; //17
    products.push(temp);
    temp = `${fminus}${minus}${bminus}`; //20
    products.push(temp);
    temp += `${plus}${fplus}`; //10
    products.push(temp);
    temp = `${bplus}${plus}${fminus}${minus}${bminus}`; //16
    products.push(temp);
    temp = `${bplus}${plus}${fplus}${minus}${bplus}${plus}${fminus}${minus}${bminus}`; // 18
    products.push(temp);
    temp += `${plus}${fplus}${minus}${bminus}${plus}${fminus}${minus}${bminus}`; //19
    products.push(temp);

    return products;
};


const enzymeCutSiteRegex = ({cutsite}) => {
    var temp = cutsite.toLowerCase();
    var reg = [];

    for (var i = 0; i < temp.length; i++){
        var str = '';

        switch (temp.charAt(i)){
            case 'a':
            case 't':
            case 'c':
            case 'g':
                str = temp.charAt(i);
                break;
            case 'm':
                str = '[ac]'; break;
            case 'k':
                str = '[gt]'; break;
            case 'r':
                str = '[ag]'; break;
            case 'y':
                str = '[ct]'; break;
            case 'w':
                str = '[at]'; break;
            case 's':
                str = '[gc]'; break;
            case 'h':
                str = '[act]'; break;
            case 'v':
                str = '[acg]'; break;
            case 'n':
                str = '[atcgmkrywshvnbd]'; break;
            case 'b':
                str = '[cgt]'; break;
            case 'd':
                str = '[agt]'; break;
            default:
                throw new Error(`error, unrecognized character ${temp.charAt(i)}`);
        }
        reg.push(str);
    }
    return new RegExp(reg.join(''), 'gi');
};


const cutProductsWithEnzyme = (products, enzyme) => {
    var re = enzymeCutSiteRegex(enzyme);
    var fragments = [];

    for (var i = 0; i < products.length; i++){
        var temp = products[i];
        var result;
        var indicies = [];

        while ((result = re.exec(temp)) !== null){
            var {index} = result;
            indicies.push(index);
        }
        var prev = 0;

        if (indicies.length == 0){
            fragments.push(products[i].length);
        }

        for (var j = 0; j < indicies.length; j++){
            var f = indicies[j] + enzyme.cutindex - prev;
            fragments.push(f);
            prev = indicies[j];
        }
    }

    if (fragments.length == 0){
        return null;
    }
    return fragments.sort();
};


const reverseCompEnzyme = (ren) => {
    var site = complementSeq(ren.cutsite, true);
    var halfway = Math.ceil(ren.cutsite.length / 2.0);

    if (site.localeCompare(ren.cutsite) && ren.cutindex == halfway){
        return null;
    }
    var result = {name: ren.name, cutindex: ren.cutsite.length - ren.cutindex, cutsite: site};
    return result;
};


/**
 *
 * @param {Array.<string>} products the LAMP product sequences
 * @param {Array.<string>} layout the lanes. A list of restriction enzyme names and ladder names
 */
const digestProducts = (products, layout) => {
    const undigestedProduct = products.map(p => p.length).sort();
    var lanes = [];

    for (const name of layout){

        if (!name) {
            lanes.push(undigestedProduct);
        } else if (name === '100bp') {
            // concat to make the ladder darker
            lanes.push(bpLadders['100'].concat(bpLadders['100']));
        } else if (name === '50bp') {
            // concat to make the ladder darker
            lanes.push(bpLadders['50'].concat(bpLadders['50']));
        } else if (enzymes[name] === undefined) {
            throw new Error(`Did not recognize restriction digest enzyme ${name}`);
        } else {
            const enzyme = enzymes[name];
            const reverseEnzyme = reverseCompEnzyme(enzyme);
            let fragments;

            if (!reverseEnzyme){
                fragments = cutProductsWithEnzyme(products, enzyme);
            } else {
                fragments = cutProductsWithEnzyme(products, enzyme);
                fragments.concat(cutProductsWithEnzyme(products, reverseEnzyme));
            }
            // if the enzyme didn't cute anything, return the undigested product
            if (fragments.length === 0){
                lanes.push(undigestedProduct);
            } else {
                lanes.push(fragments);
            }
        }
    }
    return lanes;
};


const createLampReaction = (fastaContent, lanesLayout = ['100bp', '50bp', '']) => {
    const [reaction, primerPositions] = getPrimerLayout(parseLampFasta(fastaContent));
    const blocks = createAmplificationBlocks(reaction, primerPositions);
    const products = amplifyBlocks(blocks);
    const lanes = digestProducts(products, lanesLayout);
    return {reaction, lanes, undigestedProduct: products.map(p => p.length).sort()};
};

module.exports = {
    parseLampFasta,
    createLampReaction,
    createAmplificationBlocks,
    amplifyBlocks,
    getPrimerLayout,
    complementSeq
};