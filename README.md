# demultiplexing
diagnostic existing demultiplexing issue in metagenomic virus

Command example: demultiplex match -m1 -p Result ../script/Run13tags14bases.txt *R1* *R2*

demultiplex match -mismatch_allowed -p Result_folder taglist_file fastq_R1_file fastq_R2_file

Code adapted from https://github.com/jfjlaros/demultiplex to better fit the need for viral metagnome pair-end sequencing.

