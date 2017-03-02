<template lang="pug">

div
  h1  {{ infos.title }}
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.center.
    Submit parts and a receptor vector. Get an annotated Genbank of the
    resulting construct(s).

  .form
    h4.formlabel Select an enzyme
    .enzymes-radio
      el-row(:gutter='40')
        el-col(v-for='enzyme in enzymes', :md='8', :sm='8', :xs='24')
          el-radio(v-model='form.enzyme', class="radio", :label='enzyme') {{enzyme}}
    //- h4.formlabel Provide a receptor vector
    //- sequencesuploader(v-model='form.backbone', :multiple='false',
    //-                   text="Drop a single Genbank/Fasta file (or click to select)")
    h4.formlabel Provide Parts and a receptor vector
    sequencesuploader(v-model='form.parts', :multiple='true',
                      text="Drop multiple Genbank/Fasta (or click to select)")
    results-area(:form='form', :backendUrl='infos.backendUrl')
</template>

<script>
import learnmore from '../../components/widgets/LearnMore'
import sequencesuploader from '../../components/widgets/SequencesUploader'
import digestionset from '../../components/widgets/DigestionSelectorSet'

var infos = {
  title: 'Simulate Golden Gate assemblies',
  navbarTitle: 'Simulate Golden Gate',
  path: 'simulate_gg_assemblies',
  description: '',
  backendUrl: 'start/simulate_cloning',
  icon: require('assets/images/assembly.svg')
}

export default {
  data: function () {
    return {
      enzymes: ['BsaI', 'BsmBI', 'BbsI'],
      form: {
        enzyme: 'BsaI',
        parts: [],
        backbone: null
      },
      infos: infos,
      ladder_options: [
        {
          label: 'Ladder 100 bp - 4000 bp',
          value: '100-4k'
        }
      ]
    }
  },
  components: {
    sequencesuploader,
    learnmore,
    digestionset
  },
  infos: infos,
  methods: {
    handleSuccess: function (evt) {
      console.log(evt)
    }
  }
}
</script>

<style lang='scss' scoped>

h4.formlabel {
  text-align: center;
  text-transform: uppercase;
  margin-top: 40px
}

.form {
  margin: 50px auto;
  max-width: 500px;
}

.title-img {
  height:80px;
  margin-top: -20px;
  margin-bottom: 20px;

}

.el-checkbox {
  font-weight: normal;
}

.el-select {
  width: 100%
}

.enzymes-radio {
  width: 400px;
  margin: 0 auto;
  .el-radio {
    margin-bottom: 20px;
  }
}
</style>
