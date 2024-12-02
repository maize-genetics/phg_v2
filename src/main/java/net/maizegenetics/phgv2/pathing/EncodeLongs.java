package net.maizegenetics.phgv2.pathing;

public class EncodeLongs {
    public void encodeLongs(String seq) {
        long posStrand = 0l;
        long negStrand = 0l;

        for(int i = 0; i < 31; i++) {
            char c = seq.charAt(i);
            posStrand = posStrand << 2;
            negStrand = negStrand >> 2;
            switch(c) {
                case 'A':
                    posStrand = posStrand | 0;
                    negStrand = negStrand | (3L << 62);
                    break;
                case 'C':
                    posStrand = posStrand | 1;
                    negStrand = negStrand | (2L << 62);
                    break;
                case 'G':
                    posStrand = posStrand | 2;
                    negStrand = negStrand | (1L << 62);
                    break;
                case 'T':
                    posStrand = posStrand | 3;
                    negStrand = negStrand | (0L << 62);
                    break;
                default:
                    posStrand = posStrand | 0;
                    negStrand = negStrand | (0L << 62);
                    break;
            }
        }
        for( int i = 31; i < seq.length(); i++) {
            char c = seq.charAt(i);
            posStrand = posStrand << 2;
            negStrand = negStrand >> 2;
            switch(c) {
                case 'A':
                    posStrand = posStrand | 0;
                    negStrand = negStrand | (3L << 62);
                    break;
                case 'C':
                    posStrand = posStrand | 1;
                    negStrand = negStrand | (2L << 62);
                    break;
                case 'G':
                    posStrand = posStrand | 2;
                    negStrand = negStrand | (1L << 62);
                    break;
                case 'T':
                    posStrand = posStrand | 3;
                    negStrand = negStrand | (0L << 62);
                    break;
                default:
                    posStrand = posStrand | 0;
                    negStrand = negStrand | (0L << 62);
                    break;
            }
            long minKmer = Math.min(posStrand, negStrand);
        }
    }
}
