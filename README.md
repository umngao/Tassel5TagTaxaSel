### TASSEL source repository and docs for users and developers ###

* [General overview of docs](https://bitbucket.org/tasseladmin/tassel-5-source/wiki)
* [Source](https://bitbucket.org/tasseladmin/tassel-5-source/src)
* [Latest Builds](http://www.maizegenetics.net/tassel/)

This is a tassel-5-source repository based on last commit 2017-07-22  3b1bb8010539e74ade15ec0d61fc42cc4da9fa19  

[The modified / added plugin GetSelTagTaxaDistFromDBPlugin source file is located here](https://github.com/umngao/Tassel5TagTaxaSel/blob/master/src/net/maizegenetics/analysis/gbs/v2/GetSelTagTaxaDistFromDBPlugin.java)

The command to use that plugin (after building .jar file) is this:
```bash
$tasselPath -Xms10G -Xmx20G -fork1 -GetSelTagTaxaDistFromDBPlugin \
    -db $sqlite_db \
    -o ${prj_name}_SelTagTaxaDistOutput.txt \
    -tg $seltags \
    -endPlugin -runfork1 
```
where $sqlite_db points to the database file of your GBS pipeline
-tg points to the selected alien or wheat specific tags (see [alien_2ns_predict](https://github.com/umngao/alien_2ns_predict) for examples)

By running this GetSelTagTaxaDistFromDBPlugin, users can save hours or days of compute time, and terabites of storage space for large breeding datasets.







