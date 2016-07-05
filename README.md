Align two alignments using Needleman-Wunsch algorithm.

* score is number of matching non-gap characters between two alignment columns.
	.e.g. the score  for the column pair:  AAC-ERG, ABC--RG would be 4. The
	columns 0,2,5,6, having matching, non-gap characters. Cols 1 is non-matching,
	col 3 is gaps, and col 4 has a gap and an E, so these cols, 1,3, and 4, 
	contribute 0 each to the score.
* the gap penalty going into the NW algorithm here is set in subroutine 
	Needleman-Wunsch. (At present it is zero.)

	  
