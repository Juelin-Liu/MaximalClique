# amazon=amazon
# if test -f "$amazon"; then
#     echo "$amazon exists"
# else
#     curl http://snap.stanford.edu/data/bigdata/communities/com-amazon.ungraph.txt.gz --output amazon.gz
#     gzip -d amazon.gz
#     python compact.py amazon amazon_cont
#     mv amazon_cont amazon
# fi

# youtube=youtube
# if test -f "$youtube"; then
#     echo "$youtube exists"
# else
#     curl http://snap.stanford.edu/data/bigdata/communities/com-youtube.ungraph.txt.gz --output youtube.gz
#     gzip -d youtube.gz
#     python compact.py youtube youtube_cont
#     mv youtube_cont youtube
# fi



# lj=lj
# if test -f "$lj"; then
#     echo "$lj exists"
# else
#     curl http://snap.stanford.edu/data/bigdata/communities/com-lj.ungraph.txt.gz --output lj.gz
#     gzip -d lj.gz
#     python compact.py lj lj_cont
#     mv lj_cont lj
# fi



# orkut=orkut
# if test -f "$orkut"; then
#     echo "$orkut exists"
# else
#     curl http://snap.stanford.edu/data/bigdata/communities/com-orkut.ungraph.txt.gz --output orkut.gz
#     gzip -d orkut.gz
#     python compact.py orkut orkut_cont
#     mv orkut_cont orkut
# fi

python compact.py reactome reactome_cont
mv reactome_cont reactome
cd ../src/ && make reorder;
cd ../src/ && ./reorder ../data/reactome -order gro;
cd ../src/ && ./reorder ../data/reactome -order deg;
cd ../src/ && ./reorder ../data/reactome -order degdesc;

# cd ../src/ && ./reorder ../data/youtube -order gro;
# cd ../src/ && ./reorder ../data/youtube -order deg;
# cd ../src/ && ./reorder ../data/youtube -order degdesc;

# cd ../src/ && ./reorder ../data/lj -order gro;
# cd ../src/ && ./reorder ../data/lj -order deg;
# cd ../src/ && ./reorder ../data/lj -order degdesc;

# cd ../src/ && ./reorder ../data/orkut -order gro;
# cd ../src/ && ./reorder ../data/orkut -order deg;
# cd ../src/ && ./reorder ../data/orkut -order degdesc;