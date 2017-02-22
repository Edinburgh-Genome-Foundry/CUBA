<template lang="pug">

div
  h1  {{ infos.title }}
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.center Submit a backbone and parts. Get the annotated Genbanks of all combinatorial assemblies.

  .form
    h4.formlabel Provide a backbone
    sequencesuploader(v-model='form.files', text="Drop a single file (or click to select)", :multiple='false')
    el-checkbox(v-model='form.make_report') Generate report
    resultsarea(:form='form', :backendUrl='infos.backendUrl')
</template>

<script>
import learnmore from '../../components/widgets/LearnMore'
import sequencesuploader from '../../components/widgets/SequencesUploader'
import digestionset from '../../components/widgets/DigestionSelectorSet'

var infos = {
  title: 'Sculpt A Sequence',
  navbarTitle: 'Sculpt A Sequence',
  path: 'sculpt-a-sequence',
  description: '',
  backendUrl: 'start/sculpt_a_sequence',
  icon: require('assets/images/sculpt_a_sequence.svg')
}

export default {
  data: function () {
    return {
      form: {
        ladder: '100-4k',
        digestions: [],
        make_report: false,
        files: []
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

<style scoped>

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
</style>
