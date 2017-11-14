% Inclass assignment 9
%Walter Frank Lenoir
%GB comments
1) 100
2) 100
3) 100 
4) 100
Overall: 100

% The accession number for human NOTCH1 mRNA is AF308602
% 1. Read the information from this entry into matlab
genbankinfo = getgenbank('AF308602');
% 2. Write code that runs a blast query on the first 500 base pairs of this
% gene against the refseq_rna database
seq = genbankinfo.Sequence;
seq = seq(1:500);
[requestID, requestTime] = blastncbi(seq,'blastn','Database','refseq_rna');
blast_data = getblast(requestID,'WaitTime',requestTime);
% 3. Find the three highest scoring hits from other species and identify
% the length of the alignment and fraction of matches/mismatches. 
human = blast_data.Hits(1);
% The first hit comes from Homo Sapiens NOTCH1 (expected)
chimp = blast_data.Hits(2);
lengthalign1 = length(chimp.HSPs(1).Alignment);
frac1num = chimp.HSPs(1).Identities.Match;
frac1den = chimp.HSPs(1).Identities.Possible;
% The second set of hit comes from Pan troglodytes (common chimp). It has an
% alignment length of 500 and 497 matches/ 500 possibilities.
rhino = blast_data.Hits(6);
lengthalign1 = length(rhino.HSPs(1).Alignment);
frac1num = rhino.HSPs(1).Identities.Match;
frac1den = rhino.HSPs(1).Identities.Possible;
% The third hit comes from the Rhinopithecus bieti, the black
% snub-nosed monkey. It has an alignment length of 500, and 488 matches/
% 500 possibilities. 
oldworldmonkey = blast_data.Hits(7);
lengthalign1 = length(oldworldmonkey.HSPs(1).Alignment);
frac1num = oldworldmonkey.HSPs(1).Identities.Match;
frac1den = oldworldmonkey.HSPs(1).Identities.Possible;

%The 4th hit comes from Cercocebus atys, also know as the Sooty mangabey.
%It has an alignment length of 500, and 486 matches/500 possibilities. 

% 4. Run the same query against the database est_human. Comment on the
% sequences that you find. 

[requestID, requestTime] = blastncbi(seq,'blastn','Database','est_human');
blast_data = getblast(requestID,'WaitTime',requestTime);

%The results from this query resulted in scores and alignments that were
%much shorter compared to the refseq_rna query. The results generated come
%from the homo sapiens expressed sequence tag database - containing
%libraries of cDNA clones. The first hit is the cDNA clone of the
%human NOTCH1 mRNA. Subsequent hits are DNA fragments that have similar
%sequence regions and high alignment scores.
