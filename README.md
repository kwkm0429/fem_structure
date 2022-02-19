# 構造解析+トポロジー最適化
### 概要
- FEMをC++で実装
- Eigenによる疎行列とOpenMPによる並列化で高速化

### 解析種類
- 線形静解析
- トポロジー最適化(密度法)

### 要素形状
- 1次四辺形要素

### 入力ファイル
- 節点ファイル: data-input/node.dat
- 要素ファイル: data-input/elem.dat
- 境界条件ファイル: data-input/bc.dat
- 解析設定ファイル: data-input/sim.prm

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