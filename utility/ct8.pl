#!/usr/bin/perl
###############################################
#   Operation cout tool
#   % ct.pl < file 
###############################################

#@operator=("\\+","-","\\*","\\/","**","sqrt","log","exp");
@operator=("\\+","-","/","\\*","\\*\\*","sqrt","log","exp");
@reserved=("!",   
           "program", 
           "subroutine",
           "return", 
           "end",
           "use", 
           "implicit",
           "real",      
           "integer",
           "intent",
           "do ", 
           "enddo",   
           "end do",
           "if",  
           "then",    
           "else",      
           "elseif", 
           "else if",
           "read",
           "write",   
           "print",
           "omp "
           );

#printf("%-100s", " ");
#@oplist = @operator;
#foreach $ope (@oplist) {
#  $ope =~ s/\\//g;
#  printf " %+2s", $ope;
#}
#print "\n";

& init_counter;

$if_doloop = 0;
while(<STDIN>) {
  if (/\S/) {                       # 空白文字以外

#   print $if_doloop;

    $if_rsv = 0;
    foreach $rsv (@reserved) {      # 予約語(宣言文を）サーチ
      if (/$rsv/i) {
        $if_rsv = 1;

        unless (/omp /i) {
          if ($rsv=~/do /i) {
            $if_doloop=1; 
          }
        }

#        print "--",$rsv;

        unless (/omp /i) {
          if (($rsv=~/enddo/i) or ($rsv=~/end do/i)) {
            $if_doloop=0; 
            & printsum($operator, $sum_ope);
            & init_counter;
          }
        }
      }
    } # foreach

    if ($if_rsv_old == 1) {         # 一行前が宣言文なら継続行の場合も宣言文とする
      if ($if_cont == 1) {
        $if_rsv = 1;
      }
    }
    $if_rsv_old = $if_rsv;

    if (/&$/) {                    # 文末が&なら継続行
      $if_cont = 1;
    } else {
      $if_cont = 0;
    }

    $i=0;
    foreach $ope (@operator) {     # カウンタの初期化
      $num_ope[$i]=0;
      $i++;
    }

    if ($if_rsv == 0) {
#   print;

      $i=0;
      foreach $ope (@operator) {   # 演算子のカウント
        $line = $_;
        $line = $line =~ s/$ope/$ope/g;
        $num_ope[$i]=$num_ope[$i]+$line;
#        $num_ope[$i]=$num_ope[$i]+s/$ope/$ope/g;
        $i++;
      }  # foreach
    } else {
      # nothing
    }

    if ($num_ope[4] > 0) {         # 乗算とべき乗演算の調整
      $num_ope[3] = $num_ope[3] - $num_ope[4]*2;
    }

    $line = $_;
    chomp($line);
    printf("%-100s", $line);       # ソースコードの出力
    if ($if_rsv == 0) {
      
      if ($if_doloop == 1) { 
        $i=0;
        foreach $ope (@operator) {
          printf("%3d", $num_ope[$i]);              # 演算数の出力
          $sum_ope[$i]=$sum_ope[$i]+$num_ope[$i];   # 演算数の合計を計算
          $i++;
        }  # foreach
      }
    }
    print "\n";
  }
}

#& printsum ($operator, $sum_ope);

sub init_counter {
  my ($i, $ope);

  $i=0;
  foreach $ope (@operator) {
    $sum_ope[$i]=0;
    $i++;
  }
}

sub printsum {
  my ($operator, $sum_ope) = @_;
  my ($i,$ope);

  printf("%-90s =========");
# printf("%-100s", " ");
  @oplist = @operator;
  foreach $ope (@oplist) {
    $ope =~ s/\\//g;
    printf " %+2s", $ope;
  }
  printf("==========");
  print "\n";

  printf("%-100s", " ");
  $i=0;
  foreach $ope (@operator) {
    printf("%3d", $sum_ope[$i]);
    $i++;
  }
  print "\n";
}
