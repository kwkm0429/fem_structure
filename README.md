# 構造解析+トポロジー最適化
### 概要
- FEMによる構造解析とトポロジー最適化をC++で実装
- ライブラリEigenによる疎行列を使用[1]
- OpenMPによる並列化（要素のマルチカラー化[2]を使用）

### 解析種類
- 線形静解析
- 線形座屈解析(予定)
- トポロジー最適化(密度法)

### 要素形状
- 2次元：1次四辺形要素
- 3次元：1次四面体要素、2次四面体要素（予定）

### 入力ファイル
- 節点ファイル: data-input/node.dat
- 要素ファイル: data-input/elem.dat
- 境界条件ファイル: data-input/bc.dat
- 解析設定ファイル: data-input/sim.prm
- トポロジー最適化設定ファイル：data-input/top.prm

### コンパイル
```
# OpenMPを使用しない場合
$ make
```
```
# OpenMPを使用する場合
$ make FLAGS=OPENMP
$ export OMP_NUM_THREADS=<num threads>
```

## 実行
```
$ ./bin/fem
```

### 参考文献
1. Eigen, https://eigen.tuxfamily.org/index.php?title=Main_Page
2. 有限要素計算における全体剛性行列の作成法― 疎行列データ構造の視点から ―, 永井学志, 橋本一輝, https://www.jsces.org/activity/journal/files/tutorial_1904.pdf