folders=(Demonstrators ExaHyPE Toolkit)
files=(README.md LICENSE.txt)

ExaHyPEVersion="ExaHyPE-$(git log --format="%h" -n 1 .)"

tarName="$ExaHyPEVersion.tar.gz"

echo $tarName

tar --exclude-vcs --exclude=*.o -czvf  $tarName ${folders[*]} ${files[*]}

