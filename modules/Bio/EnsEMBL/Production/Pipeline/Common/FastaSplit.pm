=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2025] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.


=head1 NAME

 Bio::EnsEMBL::Production::Pipeline::Common::FastaSplit;

=head1 DESCRIPTION

 This modules split a single fasta file into multiple chunks

=head1 MAINTAINER/AUTHOR

 jallen@ebi.ac.uk, ckong@ebi.ac.uk

=cut
package Bio::EnsEMBL::Production::Pipeline::Common::FastaSplit; 

use strict;
use warnings;
use Bio::SeqIO;
use File::Basename qw(dirname fileparse);
use File::Path qw(make_path remove_tree);
use POSIX qw(ceil);
use base ('Bio::EnsEMBL::Production::Pipeline::Common::Base');

sub param_defaults {
    my ($self) = @_;
  
    return {
      'out_dir'                 => undef,  # Top-level directory for output
      'max_seqs_per_file'       => 100,    # Maximum number of records in a file
      'max_seq_length_per_file' => undef,  # Maximum sequence length in a file
      'max_files_per_directory' => 100,    # Maximum number of files in a directory
      'max_dirs_per_directory'  => 100,    # Maximum number of subdirectories
      'delete_existing_files'   => 1,      # Ensure that directories only contains files generated with the latest execution
      'unique_file_names'       => 0,      # Ensure that output file names are unique across all directories?
      'delete_original_file'    => 0,      # After splitting, delete original fasta file?
      'file_varname'            => 'split_file',
   };
}

sub fetch_input {
    my ($self) = @_;
    my $fasta_file = $self->param_required('fasta_file');
    my $out_dir = $self->param('out_dir');
  
    if (!-e $fasta_file) {
      $self->throw("Fasta file '$fasta_file' does not exist");
    }
  
    if (!defined $out_dir) {
      $out_dir = dirname($fasta_file);
      $self->param('out_dir', $out_dir)
    } else {
     if (!-e $out_dir) {
       $self->warning("Output directory '$out_dir' does not exist. I shall create it.");
       make_path($out_dir) or $self->throw("Failed to create output directory '$out_dir'");
     }
   }
  
}

sub run {
    my $self = shift @_;
  
    # All these should have appropriate values, either from param_defaults
    # or fetch_input. (If the user has provided anything but integers for
    # the max_* parameters, the script will automatically fail when it does
    # numeric comparisons with those values, so don't bother checking them.)
    my $fasta_file     = $self->param('fasta_file');
    my $max_records    = $self->param('max_seqs_per_file');
    my $max_seq_length = $self->param('max_seq_length_per_file');
    my $max_files      = $self->param('max_files_per_directory');
    my $max_dirs       = $self->param('max_dirs_per_directory');
    $self->param('split_files', []);
  
    # Tracking parameters to determine when we need a new file or directory.
    my $record_count = 0;
    my $seq_length = 0;
    my $file_count = 1;
    my $dir_index = 0;
  
    # Do nothing if there's nothing to do...
    if (-s $fasta_file == 0) { 
      $self->input_job->autoflow(0);
      return;
    }
  
    # Need to calculate required degree of subdirectory nesting.
    $self->directory_structure();
    $self->delete_existing_files() if $self->param('delete_existing_files');
  
    my ($basename, undef, undef) = fileparse($fasta_file, qr/\.[^.]*/);
    my $split_file = $self->new_filename($dir_index, $basename, $file_count);
    my $original   = Bio::SeqIO->new(-format => 'Fasta', -file => $fasta_file);
    my $split      = Bio::SeqIO->new(-format => 'Fasta', -file => ">$split_file");
  
    while (my $seq = $original->next_seq) {
      $record_count++;
      $seq_length += $seq->length;
    
      # Checking for record_count > 1 ensures that an empty file isn't created
      # if max_records <=0 or max_seq_length < the length of the sequence.
      if (
        $record_count > 1 &&
        ( (defined($max_records) && $record_count > $max_records) ||
          (defined($max_seq_length) && $seq_length > $max_seq_length)
        )
      ) {
        $record_count = 1;
        $seq_length = $seq->length;
      
        if (defined($max_files) && $file_count >= $max_files) {
          $dir_index++;
          $file_count = 1;
        } else {
          $file_count++;
        }
        $split_file = $self->new_filename($dir_index, $basename, $file_count);
        $split = Bio::SeqIO->new(-format => 'Fasta', -file => ">$split_file");
      }
    
      my $success = $split->write_seq($seq);
      $self->throw("Failed to write sequence to '$split_file'") unless $success;
    }
  
    if ($self->param('delete_original_file')) {
      unlink $fasta_file;
    }
  
}

sub write_output {
    my ($self) = @_;
    my $file_varname = $self->param_required('file_varname');
  
    foreach my $split_file (@{$self->param('split_files')}) {
      $self->dataflow_output_id({$file_varname => $split_file}, 2);
    }
}

sub new_filename {
    my ($self, $dir_index, $basename, $file_count) = @_;
    my $out_dir = $self->param('out_dir');
    my @dirs = @{$self->param('dirs')};
    my $sub_dir = "$out_dir/".$dirs[$dir_index];
  
    if (!-e $sub_dir) {
      make_path($sub_dir) or $self->throw("Failed to create output directory '$sub_dir'");
    }
  
    my $split_file;
    if ($self->param('unique_file_names')) {
      $split_file = "$sub_dir/$basename.$dir_index.$file_count.fa";
    } else {
      $split_file = "$sub_dir/$basename.$file_count.fa";
    }
    my @split_files = (@{$self->param('split_files')}, $split_file);
    $self->param('split_files', \@split_files);
  
return $split_file;
}

sub directory_structure {
    my ($self) = @_;
  
    # This function sets an arrayref paramter with directory paths;
    # which is subsequently indexed in the new_filename function by the
    # parameter that keeps track of how many directories have been seen.
    my $max_files = $self->param('max_files_per_directory');
    my $max_dirs  = $self->param('max_dirs_per_directory');

    my $files_required = $self->files_required();
    my $dirs_required = 1;
    if (defined $max_files && $max_files > 0) {
      $dirs_required = ceil($files_required / $max_files);
    }
    if (!defined $max_dirs || $max_dirs == 0) {
      $max_dirs = 1;
    }
  
    my @dirs;
    if ($dirs_required < $max_dirs) {
      @dirs = (1..$dirs_required);
    } else {
      @dirs = (1..$max_dirs);
    }
  
   while ($dirs_required > $max_dirs) {
     $dirs_required = ceil($dirs_required / $max_dirs);
     my @new_dirs;
     foreach my $dir (@dirs) {
       foreach my $sub_dir (1..$max_dirs) {
         push @new_dirs, "$dir/$sub_dir";
       }
     }
     @dirs = @new_dirs;
  }
  
  $self->param('dirs', \@dirs);  
}

sub files_required {
    my ($self) = @_;
  
    my $fasta_file     = $self->param('fasta_file');
    my $max_records    = $self->param('max_seqs_per_file');
    my $max_seq_length = $self->param('max_seq_length_per_file');
  
    my $record_count   = 0;
    my $seq_length     = 0;
    my $files_required = 1;
  
    my $original = Bio::SeqIO->new(-format => 'Fasta', -file => $fasta_file);
  
    while (my $seq = $original->next_seq) {
      $record_count++;
      $seq_length += $seq->length;
    
      if (
        $record_count > 1 &&
        ( (defined($max_records) && $record_count > $max_records) ||
          (defined($max_seq_length) && $seq_length > $max_seq_length)
        )
      ) {
        $record_count = 1;
        $seq_length = $seq->length;
        $files_required++;
      }
   }
  
return $files_required;
}

sub delete_existing_files {
    my ($self) = @_;
  
    my $out_dir = $self->param('out_dir');
    foreach my $dir (@{$self->param('dirs')}) {
      remove_tree("$out_dir/$dir", {keep_root => 1});
    }
}

1;

