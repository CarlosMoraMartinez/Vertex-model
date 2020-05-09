


grep ">" $1 > A.txt
grep ">" $2 > B.txt
diff A.txt B.txt
rm A.txt B.txt
