cat $1 | grep -Po "(?<={\"rows\":\s).*?(?=})" |  grep -Po "(?<=\[\").*?\",\s\".*?(?=\"])" | awk 'BEGIN{FS="\", \"";OFS="\t"} {split($2,tmp,",");print $1,$2}' | sed -e "s/,//g" > $2
# find . -type f  -name "*.html" -exec ./getinfor.sh {} {}.txt \;
