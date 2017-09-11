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
