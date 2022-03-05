# 2次元片持ち梁の静解析検証

## 入力コマンド
```
$ cp ./examples/cantilever2D_FEA/*.dat ./examples/cantilever2D_FEA/*.prm  ./data-input/
$ make
$ ./bin/fem
```

## 結果
- 先端最大変位：0.797559mm（理論値：0.80mm）
![変位](displacement.png) 
- x方向応力最大値：55.7126MPa（理論値：60.0MPa)
![x方向応力](stress.png) 

## 参考
- https://www.fem-vandv.net/m12.html