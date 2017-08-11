<template lang='pug'>
.features-editor
  p.label Features minimap
  features-minimap(:featuresData='featuresData', :window='window',
                   :sequenceLength='sequence.length', @drag='updateWindowCenter',
                   @select='selectFeature')
  p.label Sequence
  sequence-viewer(:sequence='sequence', :window='window')
  features-viewer(:featuresData='featuresData', :window='window',
                  :sequenceLength='sequence.length',
                  @select='selectFeature')
  p.label Feature editor
  .editor
    el-row.location
      el-col(:xs='7', :sm='7', :lg='7')
        p Start
        el-input-number(v-model='editorFeature.start', size='small')
      el-col(:xs='7', :sm='7', :lg='7')
        p End
        el-input-number(v-model='editorFeature.end', size='small')
      el-col(:xs='5', :sm='5', :lg='5')
        p Strand
        el-select(v-model='editorFeature.strand', placeholder='Strand', size='small')
          el-option(label='None', value='0')
          el-option(label='+', value=1)
          el-option(label='-', value=-1)
    p Feature label
    el-input(v-model='editorFeature.label', placeholder='Annotation Label')
    .buttons
      el-button(size='small', icon="plus" @click='newFeature') New
      el-button(v-if='selectedFeature' size='small', icon='edit' @click='updateFeature(selectedFeature.id)') Update
      el-button(v-if='selectedFeature' size='small', icon='delete' @click='deleteFeature(selectedFeature.id)') Delete
    //- svg(:width='ntWidth*sequence.length', height='maxLines', :viewBox="'0 0 ' + sequence.length + ' ' + maxLines",
    //-     preserveAspectRatio='xMinYMin meet')
      text(v-for='nt, i in sequence', :x='i', y='1.1') {{nt}}
</template>
<script>
import sequenceviewer from './SequenceViewer'
import featuresviewer from './FeaturesViewer'
import featuresminimap from './FeaturesMinimap'

export default {
  name: 'features-editor',
  props: {
    sequence: {default: ''},
    value: {default: () => ({})}
  },
  data () {
    return {
      lineHeight: 3,
      maxLines: 8,
      // sequence: 'ATGC'.repeat(200),
      clientWidth: 500,
      zoomFactor: 15,
      windowCenter: 0,
      editorFeature: {
        id: null,
        start: 0,
        end: 1,
        strand: 1,
        label: ''
      },
      features: Object.assign({}, this.value)
      // featuresMirror: this.features
    }
  },
  components: {
    'sequence-viewer': sequenceviewer,
    'features-viewer': featuresviewer,
    'features-minimap': featuresminimap
  },
  computed: {
    featuresData: function () {
      console.log('recomputing', this.features)
      var self = this
      return Object.keys(this.features).map(function (featureID) {
        console.log(featureID, self.editorFeature.id)
        var dict = Object.assign({
          level: self.featureLevels[featureID]
        }, self.features[featureID])
        dict['selected'] = (featureID === self.editorFeature.id)
        return dict
      })
    },
    selectedFeature: function () {
      return this.features[this.editorFeature.id]
    },
    featureLevels: function () {
      var self = this
      var featureOverlaps = {}
      var featureIDs = Object.keys(this.features)
      featureIDs.map(function (id) {
        featureOverlaps[id] = []
      })
      for (var i = 0; i < featureIDs.length; i++) {
        var id1 = featureIDs[i]
        var f1 = this.features[id1]

        for (var j = i + 1; j < featureIDs.length; j++) {
          var id2 = featureIDs[j]
          var f2 = this.features[id2]
          console.log(f1.start, f1.end, f2.start, f2.end)
          if (((f1.start <= f2.start) && (f2.start < f1.end)) ||
              ((f2.start <= f1.start) && (f2.start < f2.end))) {
            featureOverlaps[id1].push(id2)
            featureOverlaps[id2].push(id1)
          }
        }
      }
      console.log('overlaps', featureOverlaps)
      var levels = {}
      var sortedFeatureIDs = featureIDs.sort(function (id1, id2) {
        var f1 = self.features[id1]
        var f2 = self.features[id2]
        return (f2.end - f2.start) - (f1.end - f1.start)
      })
      sortedFeatureIDs.map(function (id) {
        console.log(id)
        var level = 0
        while (featureOverlaps[id].some((otherID) => (levels[otherID] === level))) {
          level++
        }
        levels[id] = level
      })
      return levels
    },
    window: function () {
      var windowSize = Math.floor(this.clientWidth / this.zoomFactor)
      windowSize = windowSize + (windowSize % 2)
      var window = {
        start: this.windowCenter - windowSize / 2,
        end: this.windowCenter + windowSize / 2
      }
      if (window.start < 0) {
        window.start = 0
        window.end = windowSize
      } else if (window.end > this.sequence.length) {
        window.end = this.sequence.length
        window.start = this.sequence.length - windowSize
      }
      return window
    }
  },
  methods: {
    handleResize: function () {
      this.clientWidth = this.$el.parentElement.clientWidth
    },
    updateWindowCenter: function (c) {
      this.windowCenter = c
    },
    selectFeature: function (featureID) {
      console.log('wwwwwww', featureID)
      this.editorFeature = Object.assign({}, this.features[featureID])
    },
    newFeature: function () {
      var newID = 0
      while (this.features[newID]) {
        newID++
      }
      newID = newID.toString()
      this.$set(this.features, newID, {
        id: newID,
        start: this.editorFeature.start,
        end: this.editorFeature.end,
        strand: this.editorFeature.strand,
        label: this.editorFeature.label
      })
      this.selectFeature(newID)
    },
    updateFeature: function (featureID) {
      console.log('update', featureID)
      this.features[featureID] = Object.assign({}, this.editorFeature)
    },
    deleteFeature: function (featureID) {
      this.$delete(this.features, featureID)
    }
  },
  mounted: function () {
    window.addEventListener('resize', this.handleResize)
    this.handleResize()
  },
  watch: {
    features: {
      deep: true,
      handler: function (value) {
        this.$emit('change', value)
      }
    }
  }
}

</script>
<!-- Add "scoped" attribute to limit CSS to this component only -->
<style lang='scss' scoped>
.features-editor {
  width: 100%;
  border-radius: 5px;
  border: 5px solid #eef;
  /*padding: 8px;*/

  p.label {
    font-size: 14px;
    margin-bottom: 0;
    font-weight: bold;
  }
  .editor {
    padding: 8px;
    background-color: #fafaff;
    .buttons {
      font-size: 16px;
      padding: 3px;
      text-align: center;
    }

    p {
      margin-top: 0;
      font-size: 13px;
      margin-bottom: 0;
      margin-top: 1em;
    }
  }
}
#sequence-vuer {
    height: 300px;
}

</style>
