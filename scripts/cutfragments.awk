# Run with variables
# cutfile -- BED output file for cuts
# fragmentfile -- BED output file for fragments
BEGIN{
  FS="\t"
  OFS="\t"
}{
  strand=$6;
  read_start=$2;
  read_end=$3;
  read_id=$1;
  flag=$7;
  tlen=$11;
  if( strand == "+" ){ 
    cut_start = read_start; 
    cut_end = read_start + 1;
  } else {
    cut_start= read_end;
    cut_end = read_end + 1;
  }
  print read_id, cut_start, cut_end, "id" , "1", strand > cutfile;
  if (tlen > 0) {
    fragment_end = read_start + tlen;
    print read_id, read_start, fragment_end, ".", "." > fragmentfile;
  }
}
