=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

 Bio::EnsEMBL::Production::Pipeline::Webdatafile::WebdataFileGeneAndTranscript;

=head1 DESCRIPTION
  Compute Gene and Transcript  step for webdatafile dumps

=cut

package Bio::EnsEMBL::Production::Pipeline::Webdatafile::WebdataFileGeneAndTranscript;

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use base qw/Bio::EnsEMBL::Production::Pipeline::Common::Base/;
use Path::Tiny qw(path);
use Array::Utils qw(intersect);
use Path::Tiny qw(path);
use Carp qw/croak/;
use JSON qw/decode_json/;
use Bio::EnsEMBL::Production::Pipeline::Webdatafile::lib::GenomeLookup;
sub param_defaults {
  my ($self) = @_;
  return {
    %{$self->SUPER::param_defaults},
  };
}


sub run {

  my ($self) = @_;
  my $species = $self->param('species');
  my $current_step = $self->param('current_step') ;
  my $output = $self->param('output_path');
  my $app_path = $self->param('app_path'); 
  my $genome_data = {
    dbname     => $self->param('dbname'),
    gca        => $self->param('gca'),
    genome_id  => $self->param('genome_id'),
    species    => $self->param('species'),
    version    => $self->param('version'),
    type       => $self->param('type'),  
    root_path  => path($self->param('root_path'))
  };
  my $dba = Bio::EnsEMBL::Registry->get_DBAdaptor( $self->param('species'), $self->param('group') );
  my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor( $self->param('species'), $self->param('group'), 'Slice' );

  my $lookup = Bio::EnsEMBL::Production::Pipeline::Webdatafile::lib::GenomeLookup->new("genome_data" => $genome_data); #"root_path"=> path("/hps/nobackup2/production/ensembl/vinay/test_webdatafile"));
  my $genome = $lookup->get_genome('1');
  my $genome_report = $genome->get_genome_report();
  my @chrs = sort {$a cmp $b } map { $_->name() } grep { $_->is_assembled } @{$genome_report};
  my $transcripts_fh = $genome->genes_transcripts_path()->child('transcripts.bed')->openw;
  my $canonicals = $self->get_canonicals($genome, $dba);
  my $mane_selects = $self->get_mane_selects($genome, $dba);
  # In-memory storage for genes, transcripts and exons during processing of a single chromosome
  my $gene_record = {};
  my $transcript_record = {};
  my $exon_record = {};
  
  while(my $chr = shift @chrs) {
     $self->process_chromosome($chr, $slice_adaptor, $genome);
  }

  close $transcripts_fh;
  

}

sub process_chromosome {
  my ($self, $chr, $slice_adaptor, $genome) = @_;
  
  my $slice = $slice_adaptor->fetch_by_region( $self->param('level'), $chr );
  my $start = 1;
  my $end = $slice->length; #replaced get_length function which was using restapi
  my $chunk = $genome->chunk_size();
  my $finished = 0; 
  while($finished == 0) {
    my $new_end = ($start + $chunk)-1;
    print STDERR "Processing $chr : $start -> $new_end\n";

    if($new_end >= $end) {
      $new_end = $end;
      $finished = 1;
    }

    #process_region($chr, $start, $new_end);
    my $data = $self->get_overlaps($chr, $start, $end, $slice_adaptor);
    #update_temporary_storage($data);
    $start = $new_end+1;
  }

 
}

sub get_overlaps {
  my ($self, $chr, $start, $end, $slice_adaptor) = @_;
  my $slice = $slice_adaptor->fetch_by_region( $self->param('level'), $chr, $start, $end );
  my @data = ();
  my $genes = $slice->get_all_Genes();
  my $transcripts = $slice->get_all_Transcripts();
  my $exons = $slice->get_all_Exons();
  for my $each_gene (@$genes){

       push(@data, {
        "gene_id" => $each_gene->stable_id(),
        "source"  => $each_gene->source(),
        "logic_name" =>  "",
        "version"    => $each_gene->version(),
        "feature_type" => "gene",
        "external_name"  => $each_gene->external_name(),
        "description"  => $each_gene->description(),
        "assembly_name"  => $self->param('assembly_default'),
        "biotype"   =>  $each_gene->biotype(),
        "end"  =>  $each_gene->end(),
        "seq_region_name" => $chr,
        "strand"  => $each_gene->strand(),
        "id"  => $each_gene->stable_id(),
        "start"  => $each_gene->start()

       });

  }
 
  for my $each_transcript (@$transcripts){
    $self->warning("------------");
    for my $key (keys %$each_transcript){
        $self->warning($key);
        $self->warning($each_transcript->{$key}());
    }
    $self->warning("ens....") ;
    exit;
=head
       push(@data, {
          "source": $each_transcript->source(),
          "logic_name": "",
          "feature_type": "transcript",
          "external_name": $each_transcript->external_name(),
          "Parent": "ENSG00000146955",
          "transcript_support_level": "2",
          "seq_region_name": "7",
          "strand": 1,
          "id": "ENST00000495590",
          "transcript_id": "ENST00000495590",
          "version": 5,
          "assembly_name": "GRCh38",
          "description": null,
          "end": 140425943,
         "biotype": "protein_coding",
         "start": 140404043
       });
=cut
  }



  #my @transcripts = @{$slice->get_all_()};
   

  


  exit;
  #my $url = "https://${rest}/overlap/region/${species}/${chr}:${start}-${end}?content-type=application/json;feature=gene;feature=transcript;feature=exon;feature=cds";
}




sub get_canonicals {

  my ($self, $genome, $dba) = @_;
  my $canonicals_hash = {};
  my $result = $dba->dbc->sql_helper()->execute( -SQL =>
    'select t.stable_id as transcript_stable_id, t.stable_id from gene g join transcript t on (g.canonical_transcript_id = t.transcript_id) join seq_region sr on (t.seq_region_id = sr.seq_region_id) join coord_system cs on (sr.coord_system_id = cs.coord_system_id) where cs.species_id = '.$genome->species_id()
  );
  
  foreach my $row (@$result) {
    $canonicals_hash->{$row->[0]} = $row->[1];
  }

  return $canonicals_hash; 
  
}

sub get_mane_selects {

  my ($self, $genome, $dba) = @_;
  my $mane_hash = {};
  my $result = $dba->dbc->sql_helper()->execute( -SQL =>
    'select t.stable_id as transcript_stable_id, 1 from transcript t join transcript_attrib ta on (t.transcript_id = ta.transcript_id) join attrib_type at on (ta.attrib_type_id = at.attrib_type_id) join seq_region sr on (t.seq_region_id = sr.seq_region_id) join coord_system cs on (sr.coord_system_id = cs.coord_system_id) where cs.species_id =' . $genome->species_id()  .' and at.code = "MANE_Select"'  );

  foreach my $row (@$result) {
    $mane_hash->{$row->[0]} = $row->[1];
  }

}





sub write_output {
  my ($self) = @_;
  my $species        = $self->param('species');
  my $group          = $self->param('group');
  my $current_step   = $self->param('current_step');
  
}



1;
