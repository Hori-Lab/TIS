#!/usr/bin/perl -w

if (@ARGV != 1) 
{
   print STDOUT "USAGE: \$SCRIPT [output]\n" ;
   exit ;
}

# .F90拡張子のファイルをすべて読み込み、
# ヘッダー行のコメントからリストを作成する。
# ヘッダーの書式：
#    +-------------------------------------------------+
#    |! [subroutine name]                              |
#    |!> @brief [...... description.... ]      &       |
#    |!>  [....continued .....end]                     |
#    |                                                 |
# 先頭行は"!"とサブルーチン名のみが許される。 ただし,間のスペースの有無は自由。
# 説明の一行目に"@brief"は必須。ただし,行末の&はあってもなくてもよい。
#
# 2010/12/27 Naoto HORI

$fileOut = $ARGV[0] ;

undef(@files) ;
@files = glob("./*.F90") ;

open(OUT,">$fileOut") ;
#          123456789012345678901234567890
print OUT "SUBROUTINE NAME            | DESCRIPTION\n" ;

FILELOOP: foreach $file (@files)
{
   open(IN,"<$file") ;

   # subroutine name
   $line = <IN> ;
   if ($line =~ /^!\s*(\S*)\s*$/) 
   {
      $name =  $1 ;
   }
   else
   {
      print STDOUT "$file dose NOT have a title\n" ;
      close(IN) ;
      next ;
   }

   # first line of description
   $line = <IN> ;
   if ($line =~ /!>\s*\@brief\s*(\S.*?\S)\s*&?+\s*$/) 
   {
      $description = $1 ;
   }
   else
   {
      print STDOUT "$file dose NOT have any description\n" ;
      close(IN) ;
      next ;
   }

   while ($line = <IN>) 
   {
      # continuing lines
      if ($line =~ /!>\s*(\S.*?\S)\s*&?+\s*$/)
      {
         $description .= ' '.$1 ;
         next ;
      }
      
      # no more line
      undef (@array) ;
      @array = split(//,$description) ;
      # output subroutine name
      printf OUT "%-26s | ",$name ; 
      $break = 65 ; # 説明文一行の長さの目安
                    # この文字数を超えたら、スペースの出現とともに改行
      $ichar = 0 ;  # 文字数をカウント
      $flg_break = 0 ; # 改行直後かどうかのフラグ
      while ($char = shift @array)
      { 
         # 改行直後はスペースを出力しない
         if ($flg_break && $char eq ' ') 
         {
            next ;
         }
         $flg_break = 0;
         ++ $ichar ;
         # 改行どきか
         if ($ichar >= $break) 
         {
            # スペースなら改行する
            if ($char eq ' ')
            {
               print OUT "\n                           | " ;
               $flg_break = 1 ;
               $ichar = 0;
               next ;
            }
         }
         print OUT $char
      }
      print OUT "\n" ;
      close(IN) ;
      next FILELOOP; #次のファイルへ
   }
}


close(OUT) ;

exit ;
