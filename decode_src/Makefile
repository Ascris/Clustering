
CC = gcc

LDFLAGS = -L/usr/local/lib -lm  -lgsl


all:	decauto distauto msm

decode: Suffix.c
	$(CC)  -g MSM.c Utils.c Sequences.c MarkovEmbed.c PrefixDecode.c PrefixCode.c SuffixTree.c Suffix.c $(LDFLAGS) -o $@

test: Test.c
	$(CC)  MSM.c Utils.c Sequences.c MarkovEmbed.c SuffixTree.c PrefixDecode.c PrefixCode.c Test.c $(LDFLAGS) -o $@


distance: ComputeDistance.c
	$(CC)  -g MSM.c Utils.c Sequences.c MarkovEmbed.c SuffixTree.c PrefixDecode.c PrefixCode.c Distance.c DistSeq.c ComputeDistance.c $(LDFLAGS) -o $@

distauto: ComputeDistanceAuto.c
	$(CC)  -g MSM.c Utils.c Sequences.c MarkovEmbed.c SuffixTree.c PrefixDecode.c PrefixCode.c Distance.c DistSeq.c ComputeDistanceAuto.c $(LDFLAGS) -o $@

decauto: ComputeDecodeAuto.c
	$(CC)  -g MSM.c Utils.c Sequences.c MarkovEmbed.c SuffixTree.c PrefixDecode.c PrefixCode.c ComputeDecodeAuto.c $(LDFLAGS) -o $@

msm: ComputeMUM.c
	$(CC)  -g Utils.c Sequences.c MarkovEmbed.c SuffixTree.c Distance.c PrefixDecode.c MSM.c PrefixCode.c DistSeq.c ComputeMUM.c $(LDFLAGS) -o $@

anchor: ComputeAnchor.c
	$(CC)  -g MSM.c Utils.c Alignment.c Sequences.c MarkovEmbed.c SuffixTree.c PrefixDecode.c PrefixCode.c Distance.c DistSeq.c ComputeAnchor.c $(LDFLAGS) -o $@

relation: ShowRelation.c
	$(CC) -g MSM.c Utils.c Alignment.c Sequences.c MarkovEmbed.c SuffixTree.c PrefixDecode.c PrefixCode.c Distance.c DistSeq.c DistAln.c ShowRelation.c $(LDFLAGS) -o $@

show: ShowEvolution.c
	$(CC) -g MSM.c Utils.c Alignment.c Sequences.c MarkovEmbed.c SuffixTree.c PrefixDecode.c PrefixCode.c Distance.c DistSeq.c DistAln.c ShowEvolution.c $(LDFLAGS) -o $@

direct: CompareDirect.c
	$(CC) -g MSM.c Utils.c Alignment.c Sequences.c MarkovEmbed.c SuffixTree.c PrefixDecode.c PrefixCode.c Distance.c DistSeq.c DistAln.c CompareDirect.c $(LDFLAGS) -o $@

compare: ShowCompare.c
	$(CC) -g MSM.c Utils.c Alignment.c Sequences.c MarkovEmbed.c SuffixTree.c PrefixDecode.c PrefixCode.c Distance.c DistSeq.c DistAln.c ShowCompare.c $(LDFLAGS) -o $@
