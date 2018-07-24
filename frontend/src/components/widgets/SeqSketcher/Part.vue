<template lang='pug'>
.sketcher-part(:style='partStyle')
  .label-and-sublabel
    textarea.label(v-model='partData.label', rows=2, placeholder='(label)',
                   :class="{'hover-only': (partData.label.length === 0)}",
                   spellcheck="false")
    textarea.sublabel.grey(v-model='partData.sublabel', rows=2,
                           :class="{'hover-only': (partData.sublabel.length === 0)}",
                           placeholder='(sub-label)' spellcheck="false")
  textarea.subscript.grey(v-model='partData.subscript', placeholder='(subscript)',
                          :class="{'hover-only': (partData.subscript.length === 0)}",
                          rows=2 spellcheck="false")
  el-popover(ref='popover1', width='260', trigger='click', placement="top-start",
             v-model="popover1Visible" z-index=1000)
    .hovered-category {{hoveredCategory}}
    .categoryChoice
      .symbol-icon(v-for='category in categories', :class="category",
                 @click='partData.category = category; popover1Visible = false'
                 @mouseleave="hoveredCategory=''" @mouseover="hoveredCategory=category")


  .controls.hover-only
    //- el-popover(width='40px', trigger='click', placement="top-start",
    //-           v-model="popover2Visible")
    //-   .color-choices
    //-     .color-choice(v-for='color in colorChoices', :style="{'background-color': color}", :key='color',
    //-                 @click='partData.bg_color = color; popover2Visible = false')
    el-button-group.buttons
      el-button(@click="$emit('delete')" title='remove' size='mini' icon='el-icon-delete' circle)
      el-button(@click='partData.reversed = !partData.reversed' title='reverse' size='mini' icon='el-icon-refresh'  circle)

      el-tooltip(:enterable='true', effect='light', v-model='popover2Visible', :manual='true', title='color')
        .color-choices(slot='content')
          .color-choice(v-for='color in colorChoices', :style="{'background-color': color}", :key='color',
                      @click='partData.bg_color = color; popover2Visible = false')
        el-button(@click='popover2Visible = !popover2Visible' size='mini' icon=' el-icon-picture'  circle)
  <el-button :class="[{reversed: partData.reversed}, partData.category, 'symbol']" v-popover:popover1></el-button>
</template>

<script>
export default {
  name: 'part',
  props: {
    value: {
      default: () => ({
        category: 'user-defined',
        symbol: null,
        bg_color: null,
        label: '',
        reversed: false,
        sublabel: '',
        subscript: ''
      })
    }
  },
  data () {
    return {
      partData: Object.assign({
        category: 'user-defined',
        bg_color: null,
        label: '',
        reversed: false,
        sublabel: '',
        subscript: ''
      }, this.value),
      hoveredCategory: '',
      popover1Visible: false,
      popover2Visible: false,
      colorChoices: [
        null,
        '#f1f1fd', // blue
        '#fdf1f1', // red
        '#f2fdf1', // green
        '#fcf1fd' // purple
      ],
      categories: [
        'ATG',
        'CDS',
        'DNA-binding-element',
        'five-prime-overhang',
        'homology-arm',
        'insulator',
        'IRES',
        'ITR',
        'LTR',
        'none',
        'origin-of-replication',
        'p2A',
        'part-linker',
        'peptide-linker',
        'promoter',
        'protein-tag',
        'recombinase-recognition-sequence',
        'restriction-enzyme-recognition-site',
        'ribosome-entry-site',
        'RNA-stability-sequence',
        'terminator',
        'three-prime-overhang',
        'user-defined',
        'UTR'
      ]
    }
  },
  computed: {
    partStyle () {
      return {
        'background-color': this.partData.bg_color,
        'width': this.maxWidth
      }
    },
    maxWidth () {
      var longest = Math.max(
        this.partData.label.length,
        this.partData.sublabel.length,
        this.partData.sublabel.length
      )
      return (longest * 0.7) + 'em'
    }
  },
  watch: {
    partData: {
      deep: true,
      handler (newval) {
        this.$emit('input', newval)
      }
    }
  }
}
</script>
<!-- Add 'scoped' attribute to limit CSS to this component only -->
<style lang='scss' scoped>

$categories: (ATG CDS DNA-binding-element five-prime-overhang homology-arm
              insulator IRES ITR LTR origin-of-replication p2A part-linker
              peptide-linker promoter protein-tag none recombinase-recognition-sequence
              restriction-enzyme-recognition-site ribosome-entry-site
              RNA-stability-sequence terminator
              three-prime-overhang user-defined UTR);

.hovered-category{
  font-family: 'Raleway';
  height: 1.5em;
  text-align: center;
  color: grey;
}

.color-choices{
  .color-choice {
    display: inline-block;
    margin: 0.1em;
    width: 0.7em;
    height: 0.8em;
    border: 0.3px solid grey;
  }
}


.categoryChoice{
  text-align: center;
  .symbol-icon {
    display: inline-block;
    width: 4em;
    height: 4em;
    margin: 0.2em;
    border: 0.3px solid #eee;
    box-shadow: 1px 1px 1px #eee;
    background-size: auto 100%;
    background-repeat: no-repeat;
    background-position: 50% 100%;
    cursor: pointer;
    @each $category in $categories {
     &.#{$category} {
       background-image: url(./symbols/#{$category}.svg);
     }
   }
  }
}


.sketcher-part {
  display: inline-block;
  vertical-align: top;
  position: relative;
  margin-right: -0.3em;
  margin-left: -0.3em;
  margin-bottom: 1em;
  min-width:8em;
  max-width:15em;
  height:16em;
  font-size: 14px;
  textarea {
    border: none !important;
    outline: 0 !important;
    border-color: Transparent;
    font-family: 'Raleway';
    font-size: 1em;
    resize: none;
  }
  .hover-only {
    visibility: hidden;
    &:focus {
      visibility: visible;
    }
  }
  &:hover {
    .hover-only {
      visibility: visible;
    }
  }
  textarea, p {
   position: relative;
   font-size: 0.85em;
   text-align: center;
   background: none;
   width: 95%;
   margin-left: 2.5%;
   &.subscript {
     top: 7.5em;
   }
   &.grey {
     color: #888;
   }
  }
  .el-button.symbol {
    border: 0;
    // display: inline;
  }
  .symbol {
    position: absolute !important;
     font-size: 1em !important;

     top: 5.5em !important;
     left: 0;
     bottom: 0 !important;
     width: 100% !important;
     line-height:.8em !important;
     cursor: pointer;
     height:5em !important;
     background-size: auto 100%;
     background-repeat: no-repeat;
     background-position: 50% 100%;
     z-index: 3;
     &.reversed {
      -webkit-transform: rotate(180deg);
     }
     @each $category in $categories {
      &.#{$category} {
        background-image: url(./symbols/#{$category}.svg);
      }
    }
  }
  .controls {
    position: absolute;
    z-index: 6;
    top: 10em;
    left: 50%;
    transform: translateX(-50%);
    font-size: 1em;
    .buttons {
      position: float;
      width: 100%;
      min-width: 100px;
    }
  }
}

</style>
