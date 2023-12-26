#!/usr/bin/perl -w

#
# @File pdbToFasta_NR.pl
# @Author akande
# @Modifed by L.Cui
# @Created 15/01/2017
#This script generates pdb and fasta files for antigen and antibody

# InputFile: complex pdb
# OutputFile: antibody.pdb, antigen.pdb, antibody.fasta, antigen.fasta -> contained in a directory

use strict;
use warnings;
use Data::Dumper;
use List::MoreUtils qw(firstidx);
use List::MoreUtils qw/ uniq /;
use feature qw/say/;


our @pdbFileList;
our (%pdb_file_info,%Fasta_chainSeq,%Fasta_chainSeq_seqres,%Fasta_chainSeq_seqres1A);
our (%fastas);
our ($antibodyfasta,$antigenfasta,$antibodypdb,$antigenpdb);
our ($hc,$lc,$ac);
our %atms_info;

our %three2one = (
      'ALA' => 'A',
      'VAL' => 'V',
      'LEU' => 'L',
      'ILE' => 'I',
      'PRO' => 'P',
      'TRP' => 'W',
      'PHE' => 'F',
      'MET' => 'M',
      'GLY' => 'G',
      'SER' => 'S',
      'THR' => 'T',
      'TYR' => 'Y',
      'CYS' => 'C',
      'ASN' => 'N',
      'GLN' => 'Q',
      'LYS' => 'K',
      'ARG' => 'R',
      'HIS' => 'H',
      'ASP' => 'D',
      'GLU' => 'E',
    );
    
sub main(){
    
    #push @pdbFileList, <./*.pdb>;
    #push @pdbFileList, <$pdbDir/*.pdb>;
    push @pdbFileList, <*.pdb>;
    

    print Dumper(\@pdbFileList);
    
    
    
    foreach my $file (@pdbFileList) {
        
        print "\n $file";
    my $pname = $file;
    $pname =~ s/\.pdb$//;

    `mkdir $pname`;
    `mv $file $pname`;
    chdir($pname) or die "$!\n";
  
    # set filenames
    $antibodyfasta = $pname."_antibody.fasta";
    $antigenfasta = $pname."_antigen.fasta";
    $antibodypdb = $pname."_antibody.pdb";
    $antigenpdb = $pname."_antigen.pdb";
        
    %pdb_file_info=();
    
    &parse_pdb_files($file); #modifies %pdb_file_info
     
    ($hc,$lc,$ac)=&identifyChains();
   
    if(($hc eq '')||($lc eq '')||($ac eq '')){
       print "\n $pname has error \n";
        print $pname." hc-> ".$hc." lc->".$lc." ac->".$ac."\n";
        next();
       }
       
    #generate fasta files
    &write_to_fasta_files();
    #generate pdb files
    &write_to_pdb_files();

   }#FOREACH PDB 
   
  }#end sub

sub write_to_fasta_files(){
  my @fastas;
  push @fastas, $antibodyfasta;
  push @fastas, $antigenfasta;
  foreach my $out_map(@fastas){
   `touch $out_map`;
    open(OUTPUT, ">>$out_map") ||  ##append
   die "Unable to open $out_map for writing: $!\n";

   if($out_map =~ /antibody/){
      #write light chians
      my @lcs = split /, /, $lc;
      print OUTPUT ">light\n";
      foreach my $ch (@lcs){
        print OUTPUT $fastas{$ch}, "\n";
      }
      #write heavy chains
      my @hcs = split /, /, $hc;
      print OUTPUT ">heavy\n";
      foreach my $ch(@hcs){
        print OUTPUT $fastas{$ch}, "\n";
      }
   }elsif($out_map =~ /antigen/){
     my @ags = split /, /, $ac;
     foreach my $ch(@ags){
        print OUTPUT $fastas{$ch}, "\n";
      }
   }

    close OUTPUT ||
   die "Unable to close $out_map: $!\n"; 

     
  select()->autoflush(1);
  }  
}

sub write_to_pdb_files(){
  my @pdbs;
  push @pdbs, $antibodypdb;
  push @pdbs, $antigenpdb;
  foreach my $out_map(@pdbs){
    `touch $out_map`;
    open(OUTPUT, ">>$out_map") ||  ##append
    die "Unable to open $out_map for writing: $!\n";

    if($out_map =~ /antibody/){
      #write light chians
      my @lcs = split /, /, $lc;
      foreach my $ch (@lcs){
        print OUTPUT $atms_info{$ch};
      }
      #write heavy chains
      my @hcs = split /, /, $hc;
      foreach my $ch(@hcs){
        print OUTPUT $atms_info{$ch};
      }
   }elsif($out_map =~ /antigen/){
     my @ags = split /, /, $ac;
     foreach my $ch(@ags){
        print OUTPUT $atms_info{$ch};
      }
   }

    close OUTPUT ||
   die "Unable to close $out_map: $!\n"; 

     
  select()->autoflush(1);
  }
}

#################################
############################################
sub identifyChains{

my $icont=$pdb_file_info{'compnd'};
 

my $light="";
my $heavy="";
my $antigen="";
my $type="";

#print Dumper (\$pdb_file_info{'compnd'});
#print "CHECK ---\n";
#print Dumper(\@$icont);
 #print "\n size =".scalar(@icont)."\n";
    for (my $i = 0; $i < @$icont; $i++) {
        
        my $line= $$icont[$i];
        my $pline='';
        my $nline='';
        my $cline=$line;

        if (substr($line, 11, 6) eq "CHAIN:"){
            
            $line=substr($line,7);
            $line =~ s/\;//g;
            $line =~ s/^\s+|\s+$//g;
            my @antch = split /:/, $line;
            $line=$antch[1];
         
            $pline=$$icont[$i-1];
            if($i+1<@$icont){
            $nline=$$icont[$i+1];
            }
         
            if(($pline =~ m/HEAVY CHAIN/) or ($nline =~ m/HEAVY CHAIN/)){
                #say "heavy-->";
                
                
               $heavy=$heavy.",".$line;
            }#heavy
            elsif(($pline =~ m/LIGHT CHAIN|LAMBDA-CHAIN|KAPPA/) or ($nline =~ m/LIGHT CHAIN|LAMBDA-CHAIN|KAPPA/)){
                #say "light-->";
                $light=$light.",".$line;
            }#light
            elsif(($pline=~m/ANTIBODY/) and ($cline =~m/CHAIN: H/)){
                    #print '\n LOOK------- $line \n';
                       $heavy=$heavy.",".$line;  
                       
                   }
            elsif(($pline=~m/ANTIBODY/) and ($cline =~m/CHAIN: L/)){
               #print '\n LOOK------- $line \n';
                    $light=$light.",".$line;
                     }
                    
                   
                else{
                 $antigen .=','.$line;
                 }
        
        }#if chain
        
 

    }#for icont


#remove comma from the beginning
$heavy =~ s/^, //g;
$light =~ s/^, //g;
$antigen =~ s/^, //g;
#say " heavy chain ".$heavy."\n";
#say "light chain ".$light."\n";
#say " antigen chain ".$antigen."\n";

 return ($heavy,$light,$antigen);
}#sub identifyChains

##################################################################
#################################################################

##Split ATOM lines either on COMPND or ATOM or MODRES to identify individual chains and components
sub parse_pdb_files{
my ($seqres,$compnd,$modres,$atoms)='';

my @farray;

my $ifile=shift;
my %atms=();

my $index;

 open(my $fh, '<:encoding(UTF-8)', $ifile)
    or die qq(Could not open file "$ifile");


while (my $line = <$fh>) {
    
    if (substr($line, 0, 6) eq 'COMPND'){
            chomp $line; 
            push (@farray, $line);
    }#COMPND
    
    elsif (substr($line, 0, 4) eq 'ATOM'){
      
        #print "\n".$line."\n";
        #ATOM      1  N   GLU G  83      16.704 -42.873  87.220  1.00101.74           N 
        
        #$atoms.=$line;##store atominfo to map to individual chain pos
        
        my $serial   =  int(substr $line, 6, 5);
        my $name     =  substr $line, 12, 4;
        my $chain_id =  substr $line, 21, 1;
       #my $res_seq  =  int(substr $line, 22, 4);#position in sequence
        my $res_seq  =  substr $line, 22, 5;     ###included icode
        my $x_coord  =  substr $line, 30, 8;
        my $y_coord  =  substr $line, 38, 8;
        my $z_coord  =  substr $line, 46, 8;
        my $res_name =  substr $line, 17, 3;
        
        $res_seq=~s/\^s//g;
        
        
         my $astr=$serial.",".$x_coord.",".$y_coord.",".$z_coord.",".$name.",".$res_name.",".$res_seq;  
         $astr=~ s/\s+//g;
         #print $astr."\n";
        $atms{$chain_id}.="\n".$astr;
        $atms_info{$chain_id} .= $line;

        }#if ^ATOM
    elsif(substr($line, 0, 3) eq 'TER'){
      my $chain_id =  substr $line, 21, 1;
      $atms_info{$chain_id} .= $line;
    }

    elsif(substr($line, 0, 6) eq 'SEQRES'){
        $seqres.= $line;
        
        
     }#if SEQRES
     elsif(substr($line, 0, 6) eq 'MODRES'){
        $modres.= $line;  
     }#MODRES
        
#print $line."\n";
 
}#while



close $fh;
select()->autoflush(1);

###process chaininfo

foreach my $ch (keys %atms){
    my $chStr='';
    print "\n ----------CHAIN:$ch--------------------\n";
   
    
    my @seq=split("\n",$atms{$ch});
    my @fastapos=();
   
    for(my $i=0;$i<@seq;$i++){
         #print "\n i=$i --> $seq[$i]\n";
        if( $seq[$i] eq ''){
            
            next;
            }
       
        my @parts=split(',',$seq[$i])   ;
        if(not $parts[6]~~ @fastapos){
            #5684,74.446,-104.096,153.578,N,GLN,2
            push @fastapos,$parts[6];
            #print "\n pushed ".$parts[6]."\n";
            
            }#if
            
          my $fpos=int(firstidx { $_ eq $parts[6]} @fastapos);
         # print "\n fpos index=".$fpos."\n";
          $parts[6]=$fpos+1; ###fasta starts from 1
          my $nstr=join(',',@parts);
          $chStr.="\n".$nstr;
        }##i
    #print "\n  pdb_file_info chain $ch $chStr\n"; 
    
    
    $pdb_file_info{$ch}=$chStr; 
     
     
     #check fasta 
    my @frow=split("\n",$chStr);
    my @seqT=();
    my @seq1LT=();
        foreach my $r (@frow){
            if($r eq ''){
                
                next;}
            my @cols=split(",",$r);
            my $can=$cols[5].":".$cols[6];
            if($can ~~ @seqT){
                
                
                }#if
                
                else{
                    
                    push(@seqT,$can);
                    push(@seq1LT,$three2one{$cols[5]});
                    }
            }#foreach r
            
    #print "\n %%%%%%%%%%%$ch->".join('',@seqT)."\n";
    $fastas{$ch}=join('',@seq1LT);      
            
            
        
    }##ch


$pdb_file_info{'compnd'}=\@farray;;
$pdb_file_info{'seqres'}=$seqres;
$pdb_file_info{'modres'}=$modres;


}# sub parse_pdb


  
main();
