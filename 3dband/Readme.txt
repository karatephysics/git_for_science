・8/31
main_brillian_igor.f90 sb_30band.f90 ham.f90
brillian zoneにおける3dバンド計算
ただし、kx,ky方向のみに展開したパターン
・9/9
check_input.f90 check_diabs.f90 sb_30band.f90
固有ベクトルの規格化が1であることを高めるプログラム
dipole_main.f90 sb_dipole.f90 sb_30band.f90
双極子遷移の計算
・9/10
sb_bandinfo.f90
dipole_main.f90内の双極子計算部分などをサブルーチンに分離したもの
この後の計算は
dipole_main.f90 sb_bandinfo.f90 sb_dipole.f90 sb_30band.f90
で実行
※siliconバンドにて実装
・9/11
プログラムのツリー構造を以下の条件に変更
1.mainプログラムは３次元バンド計算を行いつつ、吸収係数導出のためのパラメータが選択的に導出できるプログラム
2.mainプログラムは、格子定数a、双極子遷移パラメータpara、バンド幅を指定ができるようなものに変更
3.dosに関してはまだ未考慮

main.f90 sb_bandinfo.f90 sb_dipole.f90 sb_30band.f90

まずはサンプルプログラムにてGeとSiで作成
main_sample.f90 sb_sample30band.f90 sb_sampleham.f90 sb_sampleparas.f90 

・9/14
9/11の続き
main.f90 sb_bandinfo.f90 sb_dipole.f90 sb_30band.f90 
の更新
とりあえず完成、格子定数を変更するとGeとSiの原子番号さえ導入すれば値がでる
pi/aを2pi/aに変更、ここはバグでした。
