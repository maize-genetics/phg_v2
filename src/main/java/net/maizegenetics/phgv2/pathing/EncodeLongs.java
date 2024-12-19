package net.maizegenetics.phgv2.pathing;

import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.longs.LongOpenHashSet;

//public class KmerIndexData {
//    int chrAndPos;
//    long gameteBitset1;
//    long getGameteBitset2;
//
//    public KmerIndexData(int chrAndPos, long gameteBitset1, long gameteBitset2) {
//        this.chrAndPos = chrAndPos;
//        this.gameteBitset1 = gameteBitset1;
//        this.getGameteBitset2 = gameteBitset2;
//    }
//
//    public int getChrAndPos() {
//        return chrAndPos;
//    }
//    public long getGameteBitset1() {
//        return gameteBitset1;
//    }
//    public long getGameteBitset2() {
//        return getGameteBitset2;
//    }
//
//    public void setChrAndPos(int chrAndPos) {
//        this.chrAndPos = chrAndPos;
//    }
//    public void setGameteBitset1(long gameteBitset1) {
//        this.gameteBitset1 = gameteBitset1;
//    }
//    public void setGameteBitset2(long gameteBitset2) {
//        this.getGameteBitset2 = gameteBitset2;
//    }
//}

public class EncodeLongs {

    static long AValue = 0;
    static long CValue = 1;
    static long GValue = 2;
    static long TValue = 3;

    //Shifting 60 for 31 mers
    static long revCompA = TValue << 60;
    static long revCompC = GValue << 60;
    static long revCompG = CValue << 60;
    static long revCompT = AValue << 60;


    public static void encodeLongs(String seq ){//, Long2ObjectOpenHashMap<KmerIndexData> keepMap, LongOpenHashSet removeSet) {
        long posStrand = 0l;
        long negStrand = 0l;

        for(int i = 0; i < 31; i++) {
            char c = seq.charAt(i);
            posStrand = posStrand << 2;
            negStrand = negStrand >> 2;
            switch(c) {
                case 'A':
                    posStrand = posStrand | AValue;
                    negStrand = negStrand | revCompA;
                    break;
                case 'C':
                    posStrand = posStrand | CValue;
                    negStrand = negStrand | revCompC;
                    break;
                case 'G':
                    posStrand = posStrand | GValue;
                    negStrand = negStrand | revCompG;
                    break;
                case 'T':
                    posStrand = posStrand | TValue;
                    negStrand = negStrand | revCompT;
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
                    posStrand = posStrand | AValue;
                    negStrand = negStrand | revCompA;
                    break;
                case 'C':
                    posStrand = posStrand | CValue;
                    negStrand = negStrand | revCompC;
                    break;
                case 'G':
                    posStrand = posStrand | GValue;
                    negStrand = negStrand | revCompG;
                    break;
                case 'T':
                    posStrand = posStrand | TValue;
                    negStrand = negStrand | revCompT;
                    break;
                default:
                    posStrand = posStrand | 0;
                    negStrand = negStrand | (0L << 62);
                    break;
            }
            //remove first 2 bits from posStrand
            posStrand = posStrand & 0x3fffffffffffffffL;
            //remove last 2 bits from negStrand
            //We actually dont need to do this because we only add to the 3rd and 4th bits from the left from and there
            //is no need for masking
//            negStrand = negStrand & 0xfffffffffffffffcL;
            long minKmer = Math.min(posStrand, negStrand);
        }
    }
}
