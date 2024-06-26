# GoMatch
スイス方式を用いる囲碁大会のための対局相手を決定するプログラム

全員が幸せになれるような対戦相手の決定方法を、組合せ最適化問題として以下のようにモデル化した。

$$ \mbox{Minimize} \sum_{i}^{N}{\\{tD(G_i)^2 + uP(G_i)^2 + vS(G_i)^2} + wE(G_i)^2\\}$$

ただし、$`G_i`$は大会参加者の中から選んだ2名の組を指し、$`D(G_i)`$、$`P(G_i)`$、$`S(G_i)`$、$`E(G_i)`$はそれぞれ、$`G_i`$の間の棋力差、既に対局した回数、勝ち星の差、身内かどうか(身内ならば1、そうでなければ0)を示す。
t、u、v、wはそれぞれ重みである。

## 使用方法

対局のグループ名ablockの場合

ablock.member 
人数をN、i番目の名前をP_i、i番目の棋力をL_iとすると
```txt
N
P_1 L_1
...
P_i L_i
...
P_N L_N
```

```bash
$ g++ main.cpp
$ ./a.out ablock
```
